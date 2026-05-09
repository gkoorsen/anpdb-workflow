"""Pass 2: for each Pass-1 "novel" ANPDB compound, find its nearest COCONUT
neighbour by Morgan-ECFP4 Tanimoto similarity.

Compounds with max Tanimoto >= THRESHOLD have a close analogue in COCONUT
(structural novelty is weaker than the InChIKey check suggested). Compounds
below the threshold remain "truly novel" even under fuzzy matching.

Caches COCONUT fingerprints to data/coconut_fps.pkl on first run; reuses them on
subsequent runs.

Outputs:
  output/anpdb_novel_tanimoto.csv  per-compound nearest neighbour + score
  output/anpdb_truly_novel.csv     subset with max_tanimoto < THRESHOLD
  output/tanimoto_summary.txt
"""

from __future__ import annotations

import csv
import io
import pickle
import sys
import time
import zipfile
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "output"

THRESHOLD = 0.85
MORGAN_RADIUS = 2
MORGAN_BITS = 2048

COCONUT_ZIP = DATA / "coconut_csv.zip"
COCONUT_FPS_CACHE = DATA / "coconut_fps.pkl"
NOVEL_CSV = OUT / "anpdb_novel.csv"


def morgan_fp(mol: Chem.Mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_BITS)


def build_coconut_fps() -> tuple[list, list[str], list[str]]:
    """Stream COCONUT CSV from the zip, compute fingerprints, return parallel lists."""
    fps: list = []
    ids: list[str] = []
    smis: list[str] = []
    n_read = n_failed = 0
    t0 = time.time()
    with zipfile.ZipFile(COCONUT_ZIP) as zf:
        name = next(n for n in zf.namelist() if n.lower().endswith(".csv"))
        print(f"  reading {name}", file=sys.stderr)
        with zf.open(name) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(text)
            for row in reader:
                n_read += 1
                smi = (row.get("canonical_smiles") or "").strip()
                if not smi:
                    n_failed += 1
                    continue
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    n_failed += 1
                    continue
                fps.append(morgan_fp(mol))
                ids.append(row.get("identifier") or "")
                smis.append(smi)
                if n_read % 50_000 == 0:
                    rate = n_read / max(time.time() - t0, 1e-6)
                    print(
                        f"  COCONUT FPs: read={n_read:,} kept={len(fps):,} "
                        f"failed={n_failed:,} ({rate:,.0f}/s)",
                        file=sys.stderr,
                    )
    print(f"  COCONUT total kept: {len(fps):,} / read {n_read:,}", file=sys.stderr)
    return fps, ids, smis


def load_or_build_coconut() -> tuple[list, list[str], list[str]]:
    if COCONUT_FPS_CACHE.exists():
        print(f"Loading cached COCONUT fingerprints from {COCONUT_FPS_CACHE}", file=sys.stderr)
        with open(COCONUT_FPS_CACHE, "rb") as fh:
            fps, ids, smis = pickle.load(fh)
        print(f"  loaded {len(fps):,} fingerprints", file=sys.stderr)
        return fps, ids, smis
    print("Building COCONUT fingerprints (first run; will be cached) ...", file=sys.stderr)
    fps, ids, smis = build_coconut_fps()
    print(f"Caching to {COCONUT_FPS_CACHE} ...", file=sys.stderr)
    with open(COCONUT_FPS_CACHE, "wb") as fh:
        pickle.dump((fps, ids, smis), fh, protocol=pickle.HIGHEST_PROTOCOL)
    return fps, ids, smis


def main() -> int:
    print(f"Loading Pass-1 novel set from {NOVEL_CSV}", file=sys.stderr)
    novel = pd.read_csv(NOVEL_CSV, dtype=str, keep_default_na=False)
    print(f"  {len(novel):,} novel ANPDB compounds", file=sys.stderr)

    print("Computing ANPDB query fingerprints ...", file=sys.stderr)
    q_fps: list = []
    q_keep: list[int] = []
    for i, smi in enumerate(novel["smiles"]):
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None:
            q_fps.append(morgan_fp(mol))
            q_keep.append(i)
    print(f"  query FPs: {len(q_fps):,} (skipped {len(novel) - len(q_fps)} unparseable SMILES)",
          file=sys.stderr)

    coconut_fps, coconut_ids, coconut_smis = load_or_build_coconut()

    print(f"Searching {len(q_fps):,} queries against {len(coconut_fps):,} COCONUT FPs ...",
          file=sys.stderr)
    max_sim = [0.0] * len(novel)
    nearest_id = [""] * len(novel)
    nearest_smi = [""] * len(novel)
    t0 = time.time()
    for k, (orig_idx, qfp) in enumerate(zip(q_keep, q_fps)):
        sims = DataStructs.BulkTanimotoSimilarity(qfp, coconut_fps)
        # argmax without numpy
        best_i = 0
        best_s = sims[0]
        for j in range(1, len(sims)):
            s = sims[j]
            if s > best_s:
                best_s = s
                best_i = j
        max_sim[orig_idx] = best_s
        nearest_id[orig_idx] = coconut_ids[best_i]
        nearest_smi[orig_idx] = coconut_smis[best_i]
        if (k + 1) % 50 == 0:
            elapsed = time.time() - t0
            done = k + 1
            eta = elapsed / done * (len(q_fps) - done)
            print(f"  queries done: {done:,}/{len(q_fps):,}  "
                  f"elapsed={elapsed:,.0f}s  eta={eta:,.0f}s",
                  file=sys.stderr)

    novel["max_tanimoto"] = max_sim
    novel["nearest_coconut_id"] = nearest_id
    novel["nearest_coconut_smiles"] = nearest_smi
    novel["has_close_analogue"] = novel["max_tanimoto"] >= THRESHOLD
    novel["truly_novel"] = ~novel["has_close_analogue"]

    out_full = OUT / "anpdb_novel_tanimoto.csv"
    out_truly = OUT / "anpdb_truly_novel.csv"
    novel.to_csv(out_full, index=False)
    novel[novel["truly_novel"]].to_csv(out_truly, index=False)

    truly_ids_path = OUT / "truly_novel_molecule_ids.txt"
    with open(truly_ids_path, "w") as fh:
        for mid in novel.loc[novel["truly_novel"], "molecule_id"]:
            fh.write(f"{mid}\n")

    n_total = len(novel)
    n_close = int(novel["has_close_analogue"].sum())
    n_truly = int(novel["truly_novel"].sum())
    bins = [0.0, 0.50, 0.70, 0.80, 0.85, 0.90, 0.95, 1.001]
    bin_labels = ["<0.50", "0.50-0.70", "0.70-0.80", "0.80-0.85",
                  "0.85-0.90", "0.90-0.95", "0.95-1.00"]
    counts = pd.cut(novel["max_tanimoto"], bins=bins, right=False, labels=bin_labels).value_counts().sort_index()

    lines = [
        f"Pass-1 residual (no exact match in COCONUT or ChEMBL): {n_total:,}",
        f"Morgan radius={MORGAN_RADIUS}  bits={MORGAN_BITS}  threshold={THRESHOLD}",
        "",
        f"Close analogue in COCONUT (max Tanimoto >= {THRESHOLD}): {n_close:,}",
        f"Truly novel even under similarity                        : {n_truly:,}",
        "",
        "Distribution of max Tanimoto vs COCONUT:",
    ]
    for label, count in counts.items():
        lines.append(f"  {label:<12} {count:>5}")
    lines += [
        "",
        f"Outputs:",
        f"  {out_full}",
        f"  {out_truly}",
        f"  {truly_ids_path}",
    ]
    summary_path = OUT / "tanimoto_summary.txt"
    summary_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
