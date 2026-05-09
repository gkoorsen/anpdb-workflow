"""Local ChEMBL nearest-neighbour target inference for the 1,012 truly-novel
ANPDB compounds.

Method
------
For each query compound q in the novel set:
  1. compute Morgan-ECFP4 fingerprint
  2. find top-K most similar bioactive ChEMBL compounds by Tanimoto
  3. each neighbour n contributes target_score(t) += sim(q,n) for every target
     t against which n is recorded as ACTIVE at the chosen activity threshold
     in PIDGINv4's bioactivity dataset
  4. aggregate, normalise by sum-similarity, output ranked target list per
     compound

Active-evidence universe = PIDGIN's no_ortho/bioactivity_dataset (3,371 per-
target TSVs, columns include 0.1_Flag/1_Flag/10_Flag/100_Flag indicating
activity at 100 nM / 1 uM / 10 uM / 100 uM). Default threshold = 0.1 (100 nM)
to match the existing PIDGIN ad60/ad90 hits we want to compare against.

Outputs (under output/chembl_nn/)
---------------------------------
- chembl_actives_index.pkl       (built once, cached)
    {
      "compound_smiles": {chembl_id: smiles},
      "compound_targets": {chembl_id: set[uniprot]},
      "all_fps": [ExplicitBitVect, ...],
      "all_ids": [chembl_id, ...],
    }
- chembl_nn_predictions.tsv      one row per (compound, target) with score,
                                 supporting neighbour count, top neighbour ID
- chembl_nn_summary.txt
"""

from __future__ import annotations

import csv
import io
import pickle
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
csv.field_size_limit(1 << 30)

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DATA = ROOT / "data"
OUT = ROOT / "output" / "chembl_nn"
OUT.mkdir(parents=True, exist_ok=True)

PIDGIN_BIOACT = Path("/Users/gerritkoorsen/PIDGINv4/no_ortho/bioactivity_dataset")
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"

ACTIVE_THRESHOLD = "0.1"          # 100 nM, matches PIDGIN runs we want to compare
ACTIVE_FLAG_COL  = f"{ACTIVE_THRESHOLD}_Flag"
TOP_K            = 20             # number of nearest neighbours per query
SIM_FLOOR        = 0.30           # ignore neighbours below this Tanimoto
MIN_NEIGHBOUR_SUPPORT = 2         # require >= 2 supporting actives per target
MIN_NORMALISED_SCORE  = 0.40      # report rows with normalised score >= this

CACHE = DATA / "chembl_actives_index.pkl"

MORGAN_RADIUS = 2
MORGAN_BITS = 2048


def morgan_fp(mol: Chem.Mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_BITS)


def build_index() -> dict:
    """Stream every per-target zip in PIDGIN's bioactivity_dataset, collect
    unique active ChEMBL compounds at the chosen threshold, and fingerprint them.
    """
    files = sorted(PIDGIN_BIOACT.glob("*.smi.zip"))
    print(f"Scanning {len(files):,} bioactivity zips ...", file=sys.stderr)

    smiles_of: dict[str, str] = {}
    targets_of: dict[str, set[str]] = defaultdict(set)
    n_rows = 0
    n_active = 0
    t0 = time.time()

    for i, zpath in enumerate(files, start=1):
        try:
            with zipfile.ZipFile(zpath) as zf:
                inner = next(n for n in zf.namelist() if n.endswith(".smi"))
                with zf.open(inner) as raw:
                    text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
                    reader = csv.DictReader(text, delimiter="\t")
                    if ACTIVE_FLAG_COL not in (reader.fieldnames or []):
                        continue
                    for row in reader:
                        n_rows += 1
                        if row.get(ACTIVE_FLAG_COL) != "1":
                            continue
                        cid = row.get("Compound_ID", "")
                        smi = row.get("Smiles", "")
                        uni = row.get("Uniprot", "")
                        if not cid or not smi or not uni:
                            continue
                        n_active += 1
                        smiles_of.setdefault(cid, smi)
                        targets_of[cid].add(uni)
        except Exception as exc:
            print(f"  WARN: skipping {zpath.name}: {exc}", file=sys.stderr)
        if i % 200 == 0:
            elapsed = time.time() - t0
            print(f"  zips read: {i:,}/{len(files):,}  rows={n_rows:,}  "
                  f"active rows={n_active:,}  unique cpds={len(smiles_of):,}  "
                  f"({elapsed:,.0f}s)", file=sys.stderr)

    print(f"Indexed {len(smiles_of):,} unique active compounds across "
          f"{sum(len(t) for t in targets_of.values()):,} (compound, target) actives at "
          f"{ACTIVE_THRESHOLD} uM threshold", file=sys.stderr)

    print("Computing Morgan fingerprints ...", file=sys.stderr)
    all_ids: list[str] = []
    all_fps: list = []
    n_failed = 0
    t0 = time.time()
    for j, (cid, smi) in enumerate(smiles_of.items(), start=1):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_failed += 1
            continue
        all_ids.append(cid)
        all_fps.append(morgan_fp(mol))
        if j % 50_000 == 0:
            elapsed = time.time() - t0
            print(f"  fps: {j:,}/{len(smiles_of):,}  fails={n_failed}  "
                  f"({elapsed:,.0f}s)", file=sys.stderr)

    index = {
        "compound_smiles": smiles_of,
        "compound_targets": dict(targets_of),
        "all_ids": all_ids,
        "all_fps": all_fps,
        "active_threshold_uM": ACTIVE_THRESHOLD,
    }
    print(f"  cached fingerprints: {len(all_fps):,} of {len(smiles_of):,}",
          file=sys.stderr)
    return index


def get_or_build_index() -> dict:
    if CACHE.exists():
        print(f"Loading cached index from {CACHE} ...", file=sys.stderr)
        with open(CACHE, "rb") as fh:
            ix = pickle.load(fh)
        print(f"  {len(ix['all_ids']):,} active ChEMBL compounds", file=sys.stderr)
        return ix
    ix = build_index()
    print(f"Caching index to {CACHE} ({sum(1 for _ in ix['all_ids']):,} cpds) ...",
          file=sys.stderr)
    with open(CACHE, "wb") as fh:
        pickle.dump(ix, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return ix


def main() -> int:
    ix = get_or_build_index()
    chembl_ids = ix["all_ids"]
    chembl_fps = ix["all_fps"]
    cpd_targets = ix["compound_targets"]

    print(f"\nLoading novel queries from {NOVEL_STD}", file=sys.stderr)
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    print(f"  {len(novel):,} compounds", file=sys.stderr)

    # Build query fingerprints
    print("Computing query fingerprints ...", file=sys.stderr)
    q_ids: list[str] = []
    q_fps: list = []
    n_q_fail = 0
    for _, r in novel.iterrows():
        smi = r["std_smiles"] or r["smiles"]
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            n_q_fail += 1
            continue
        q_ids.append(r["molecule_id"])
        q_fps.append(morgan_fp(mol))
    print(f"  {len(q_fps):,} query FPs (skipped {n_q_fail})", file=sys.stderr)

    # k-NN target inference
    print(f"\nRunning k-NN inference (k={TOP_K}, sim_floor={SIM_FLOOR}) ...",
          file=sys.stderr)
    out_rows: list[dict] = []
    t0 = time.time()
    for i, (qid, qfp) in enumerate(zip(q_ids, q_fps)):
        sims = DataStructs.BulkTanimotoSimilarity(qfp, chembl_fps)
        # Get top-K above sim_floor
        idxs = sorted(range(len(sims)), key=lambda j: -sims[j])[:TOP_K * 4]
        kept = [(j, sims[j]) for j in idxs if sims[j] >= SIM_FLOOR][:TOP_K]
        if not kept:
            continue
        # Aggregate per target
        per_target_score: dict[str, float] = defaultdict(float)
        per_target_count: dict[str, int] = defaultdict(int)
        per_target_top: dict[str, tuple[str, float]] = {}
        sum_sim = sum(s for _, s in kept)
        for j, s in kept:
            cid = chembl_ids[j]
            for tgt in cpd_targets.get(cid, ()):
                per_target_score[tgt] += s
                per_target_count[tgt] += 1
                if tgt not in per_target_top or per_target_top[tgt][1] < s:
                    per_target_top[tgt] = (cid, s)
        for tgt, sc in per_target_score.items():
            n = per_target_count[tgt]
            top_cid, top_sim = per_target_top[tgt]
            normalised = sc / sum_sim if sum_sim > 0 else 0.0
            if n < MIN_NEIGHBOUR_SUPPORT or normalised < MIN_NORMALISED_SCORE:
                continue
            out_rows.append({
                "compound_id": qid,
                "uniprot": tgt,
                "n_supporting_neighbours": n,
                "raw_score": round(sc, 4),
                "normalised_score": round(normalised, 4),
                "top_neighbour_chembl_id": top_cid,
                "top_neighbour_tanimoto": round(top_sim, 4),
            })
        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (len(q_fps) - i - 1)
            print(f"  q={i+1:,}/{len(q_fps):,}  rows so far={len(out_rows):,}  "
                  f"elapsed={elapsed:,.0f}s  eta={eta:,.0f}s",
                  file=sys.stderr)

    pred_path = OUT / "chembl_nn_predictions.tsv"
    pd.DataFrame(out_rows).to_csv(pred_path, sep="\t", index=False)
    print(f"\nWrote {len(out_rows):,} predicted (compound, target) pairs to {pred_path}",
          file=sys.stderr)

    # Summary
    df = pd.DataFrame(out_rows)
    n_cpds = df["compound_id"].nunique() if len(df) else 0
    n_tgts = df["uniprot"].nunique() if len(df) else 0
    summary_path = OUT / "chembl_nn_summary.txt"
    with open(summary_path, "w") as fh:
        fh.write(f"ChEMBL nearest-neighbour target inference\n")
        fh.write(f"=========================================\n\n")
        fh.write(f"Active-compound universe: PIDGIN bioactivity_dataset @ "
                 f"{ACTIVE_THRESHOLD} uM threshold\n")
        fh.write(f"Active ChEMBL compounds:  {len(chembl_ids):,}\n")
        fh.write(f"Distinct target universe: "
                 f"{len({t for ts in cpd_targets.values() for t in ts}):,}\n\n")
        fh.write(f"Query compounds:           {len(q_fps):,}\n")
        fh.write(f"k-NN parameters:           K={TOP_K}, sim_floor={SIM_FLOOR}, "
                 f"min_support={MIN_NEIGHBOUR_SUPPORT}, "
                 f"min_score={MIN_NORMALISED_SCORE}\n\n")
        fh.write(f"Predictions:               {len(out_rows):,} "
                 f"({n_cpds:,} cpds, {n_tgts:,} targets)\n")
    print(open(summary_path).read())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
