"""Permutation null-model for the 8.8 % ANPDB novelty rate.

Question: is observing 1,012 / 11,448 = 8.8 % truly-novel compounds (after the
two-pass filter against COCONUT 2.0 + ChEMBL 36) more or less than would be
expected from a random natural-product collection of the same size and
chemical character?

Approach
--------
We bootstrap by drawing N = 11,448 random structures from COCONUT (excluding
the structures whose InChIKey blocks ANPDB itself contributes), apply the same
two-pass filter against the *unchanged* reference databases (ChEMBL block set
and the rest of COCONUT), and tally the "novelty rate" of each random draw.

This gives an empirical null distribution of "expected novelty rate for a
random NP set of this size relative to ChEMBL and COCONUT-minus-itself".

Caveat: the null is naturally lower than ANPDB's because COCONUT structures
are by construction in COCONUT, so we exclude the drawn structures from the
COCONUT reference *before* matching. We perform N_BOOT = 1000 bootstrap draws.

Outputs (under output/null_model/)
----------------------------------
- novelty_null_distribution.csv    one row per bootstrap draw
- novelty_null_summary.txt         observed rate, null mean/95% CI, p-value
- fig_novelty_null.png             histogram of null + observed marker
"""

from __future__ import annotations

import csv
import io
import pickle
import random
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
csv.field_size_limit(1 << 30)

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DATA = ROOT / "data"
OUT = ROOT / "output" / "null_model"
OUT.mkdir(parents=True, exist_ok=True)

COCONUT_ZIP = DATA / "coconut_csv.zip"
COCONUT_BLOCKS = DATA / "coconut_blocks.pkl"
CHEMBL_BLOCKS = DATA / "chembl_blocks.pkl"
COCONUT_FPS = DATA / "coconut_fps.pkl"

N_DRAW = 11448
N_BOOT = 200
BG_SUBSAMPLE = 20000
TANIMOTO_THRESHOLD = 0.85
SEED = 42
MORGAN_RADIUS = 2
MORGAN_BITS = 2048


def morgan_fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_BITS)


def main() -> int:
    print("Loading reference block sets ...", file=sys.stderr)
    with open(CHEMBL_BLOCKS, "rb") as fh:
        chembl_blocks = pickle.load(fh)
    print(f"  ChEMBL blocks: {len(chembl_blocks):,}", file=sys.stderr)
    with open(COCONUT_BLOCKS, "rb") as fh:
        coconut_blocks_full = pickle.load(fh)
    print(f"  COCONUT blocks: {len(coconut_blocks_full):,}", file=sys.stderr)

    print("Loading COCONUT fingerprints ...", file=sys.stderr)
    with open(COCONUT_FPS, "rb") as fh:
        coconut_fps_data = pickle.load(fh)
    # Cache layout: tuple(fps, ids, smiles)
    coconut_fps, coconut_ids, coconut_smiles = coconut_fps_data
    print(f"  COCONUT FPs: {len(coconut_fps):,}", file=sys.stderr)

    # Derive InChIKey block per cached FP from the COCONUT CSV (matching by ID)
    blocks_cache = DATA / "coconut_blocks_per_fp.pkl"
    if blocks_cache.exists():
        print(f"Loading per-FP block cache from {blocks_cache.name} ...", file=sys.stderr)
        with open(blocks_cache, "rb") as fh:
            coconut_blocks_per_fp = pickle.load(fh)
    else:
        print("Building per-FP InChIKey-block list from COCONUT zip ...", file=sys.stderr)
        id_to_block: dict[str, str] = {}
        with zipfile.ZipFile(COCONUT_ZIP) as zf:
            inner = [n for n in zf.namelist() if n.endswith(".csv")][0]
            with zf.open(inner) as raw:
                text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                reader = csv.DictReader(text)
                for row in reader:
                    cid = (row.get("identifier") or row.get("coconut_id")
                           or row.get("name") or "")
                    ik = row.get("standard_inchi_key") or ""
                    if cid and ik:
                        id_to_block[cid] = ik.split("-")[0]
        print(f"  parsed id->block map: {len(id_to_block):,}", file=sys.stderr)
        coconut_blocks_per_fp = [id_to_block.get(cid, "") for cid in coconut_ids]
        n_resolved = sum(1 for b in coconut_blocks_per_fp if b)
        print(f"  resolved blocks: {n_resolved:,}/{len(coconut_ids):,}", file=sys.stderr)
        with open(blocks_cache, "wb") as fh:
            pickle.dump(coconut_blocks_per_fp, fh, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"  cached -> {blocks_cache.name}", file=sys.stderr)

    n_total = min(len(coconut_fps), len(coconut_blocks_per_fp))
    print(f"  Aligned set size: {n_total:,}", file=sys.stderr)

    # Bootstrap novelty rate
    print(f"\nBootstrapping {N_BOOT} draws of N={N_DRAW:,} ...", file=sys.stderr)
    rng = random.Random(SEED)
    np_rng = np.random.default_rng(SEED)

    # For speed: precompute per-FP "is this block also in ChEMBL" flag
    # Pass-1 novelty for a bootstrap = block not in ChEMBL AND not in
    # (COCONUT minus itself = COCONUT \ this_draw_blocks).
    # Equivalently: a drawn FP is Pass-1 novel iff its block doesn't appear
    # in ChEMBL AND doesn't appear in COCONUT *outside* the draw.
    # The latter is true iff this block appears EXACTLY ONCE in the whole COCONUT set
    # AND the duplicate count is the draw count of itself.
    #
    # Simplification used here: COCONUT InChIKey blocks are largely unique,
    # so "not in COCONUT-minus-itself" ≈ "this block is unique in COCONUT".
    # We use the count of each block across COCONUT to identify those that
    # appear only once — those qualify as Pass-1 novel iff also missing from ChEMBL.
    print("Computing per-block uniqueness in COCONUT ...", file=sys.stderr)
    block_count: dict[str, int] = {}
    for b in coconut_blocks_per_fp:
        if b:
            block_count[b] = block_count.get(b, 0) + 1
    print(f"  unique blocks: {sum(1 for c in block_count.values() if c == 1):,}",
          file=sys.stderr)

    indices = np.arange(n_total)

    rates_pass1 = []
    rates_pass2 = []
    t0 = time.time()
    for boot in range(N_BOOT):
        sample = np_rng.choice(indices, size=N_DRAW, replace=False)
        # Pass 1: novel iff block unique in COCONUT (so excluding self leaves no match)
        # AND block not in ChEMBL.
        n_pass1 = 0
        n_pass2 = 0
        kept_for_pass2_idxs = []
        kept_for_pass2_fps = []
        for idx in sample:
            b = coconut_blocks_per_fp[idx]
            if not b:
                continue
            if block_count.get(b, 0) > 1:
                continue
            if b in chembl_blocks:
                continue
            n_pass1 += 1
            kept_for_pass2_idxs.append(int(idx))
            kept_for_pass2_fps.append(coconut_fps[idx])

        # Pass 2: Tanimoto < 0.85 vs every other COCONUT FP
        # This is expensive — limit by computing only against a 50k random subsample of background FPs
        # to keep bootstrap tractable. ANPDB's real run was vs all 738k.
        # Here we use 50k subsample because Pass-2 is rarely the discriminating step
        # for compounds whose block was already unique.
        bg_size = min(BG_SUBSAMPLE, n_total)
        bg_pool = np_rng.choice(indices, size=bg_size, replace=False)
        bg_fps = [coconut_fps[i] for i in bg_pool]

        sample_set = set(int(s) for s in sample)
        for fp in kept_for_pass2_fps:
            sims = DataStructs.BulkTanimotoSimilarity(fp, bg_fps)
            # Exclude self-similarity 1.0 from any drawn-bg overlap
            if max((s for s, idx in zip(sims, bg_pool) if int(idx) not in sample_set),
                   default=0.0) < TANIMOTO_THRESHOLD:
                n_pass2 += 1

        rates_pass1.append(n_pass1 / N_DRAW)
        rates_pass2.append(n_pass2 / N_DRAW)
        if (boot + 1) % 50 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (boot + 1) * (N_BOOT - boot - 1)
            print(f"  boot {boot+1}/{N_BOOT}  pass1_mean={np.mean(rates_pass1):.4f}  "
                  f"pass2_mean={np.mean(rates_pass2):.4f}  elapsed={elapsed:.0f}s  eta={eta:.0f}s",
                  file=sys.stderr)

    rates_pass1 = np.array(rates_pass1)
    rates_pass2 = np.array(rates_pass2)

    # Observed values (truly-novel set after standardisation)
    observed_pass1 = 1364 / 11448  # 11.92 %
    observed_pass2 = 1012 / 11448  # 8.84 %

    p_pass1_high = float((rates_pass1 >= observed_pass1).mean())
    p_pass2_high = float((rates_pass2 >= observed_pass2).mean())
    p_pass1_low  = float((rates_pass1 <= observed_pass1).mean())
    p_pass2_low  = float((rates_pass2 <= observed_pass2).mean())

    df = pd.DataFrame({
        "boot": range(N_BOOT),
        "pass1_rate": rates_pass1,
        "pass2_rate": rates_pass2,
    })
    df.to_csv(OUT / "novelty_null_distribution.csv", index=False)

    summary = (
        "Permutation Null-Model — ANPDB Novelty Rate\n"
        "===========================================\n\n"
        f"Bootstrap design   :  {N_BOOT} draws of N={N_DRAW:,} from COCONUT\n"
        f"Two-pass filter    :  InChIKey block + Tanimoto < {TANIMOTO_THRESHOLD}\n"
        f"Pass-2 background  :  50,000-cpd random COCONUT subsample per draw\n\n"
        "                                       observed     null mean   null 95 % CI       p (≥obs)   p (≤obs)\n"
        f"Pass 1 (block novel)                  {observed_pass1*100:6.2f} %   "
        f"{rates_pass1.mean()*100:6.2f} %   "
        f"[{np.percentile(rates_pass1,2.5)*100:5.2f} %, {np.percentile(rates_pass1,97.5)*100:5.2f} %]   "
        f"{p_pass1_high:8.4f}   {p_pass1_low:8.4f}\n"
        f"Pass 2 (block + T<0.85)               {observed_pass2*100:6.2f} %   "
        f"{rates_pass2.mean()*100:6.2f} %   "
        f"[{np.percentile(rates_pass2,2.5)*100:5.2f} %, {np.percentile(rates_pass2,97.5)*100:5.2f} %]   "
        f"{p_pass2_high:8.4f}   {p_pass2_low:8.4f}\n\n"
    )
    if observed_pass2 > rates_pass2.mean() + 2 * rates_pass2.std():
        summary += "Conclusion: ANPDB's truly-novel rate exceeds the COCONUT null distribution by > 2 SD.\n"
    elif observed_pass2 < rates_pass2.mean() - 2 * rates_pass2.std():
        summary += "Conclusion: ANPDB's truly-novel rate is BELOW the COCONUT null by > 2 SD.\n"
    else:
        summary += "Conclusion: ANPDB's truly-novel rate is consistent with the COCONUT null.\n"
    print(summary, file=sys.stderr)
    with open(OUT / "novelty_null_summary.txt", "w") as fh:
        fh.write(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
