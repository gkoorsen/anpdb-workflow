"""Tightened novelty: standardize ANPDB structures (largest fragment, normalize,
uncharge), recompute parent InChIKey-skeleton blocks, look them up in COCONUT
and ChEMBL, then run Pass 2 (Morgan-ECFP4 Tanimoto vs COCONUT) on the residual.

COCONUT 2.0 is already standardized via the ChEMBL Curation Pipeline, so
applying standardization only to the ANPDB side rescues "ANPDB salt vs COCONUT
parent" duplicates.

Reference InChIKey-block sets and COCONUT fingerprints are loaded from the
caches written by earlier passes; rebuild them if missing.

Outputs:
  output/anpdb_annotated_std.csv          every ANPDB row + standardized SMILES,
                                          parent InChIKey block, and in_*/novel
                                          flags for raw and standardized passes
  output/anpdb_novel_std.csv              Pass-1 standardized novel residual
  output/anpdb_truly_novel_std.csv        Pass-2 (Tanimoto) standardized novel
  output/novel_molecule_ids_std.txt       Pass-1 standardized IDs
  output/truly_novel_molecule_ids_std.txt Pass-2 standardized IDs (TIGHTENED)
  output/summary_std.txt                  before/after counts
"""

from __future__ import annotations

import csv
import gzip
import io
import pickle
import re
import sys
import time
import zipfile
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "output"

ANPDB_PATH = DATA / "ANPDB.csv"
COCONUT_ZIP = DATA / "coconut_csv.zip"
CHEMBL_PATH = DATA / "chembl_36_chemreps.txt.gz"
COCONUT_BLOCKS_CACHE = DATA / "coconut_blocks.pkl"
CHEMBL_BLOCKS_CACHE = DATA / "chembl_blocks.pkl"
COCONUT_FPS_CACHE = DATA / "coconut_fps.pkl"

THRESHOLD = 0.85
MORGAN_RADIUS = 2
MORGAN_BITS = 2048

_normalizer = rdMolStandardize.Normalizer()
_largest = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()
_sep_re = re.compile(r" {4,}")


def standardize(smi: str | None, inchi: str | None = None) -> Chem.Mol | None:
    """Strip salts, normalize, uncharge. Tries SMILES first, then InChI as
    fallback (for ANPDB rows missing SMILES).
    """
    mol = None
    if smi:
        mol = Chem.MolFromSmiles(smi)
    if mol is None and inchi:
        mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return None
    try:
        mol = _normalizer.normalize(mol)
        mol = _largest.choose(mol)
        mol = _uncharger.uncharge(mol)
    except Exception:
        return None
    return mol


def inchikey_block(mol: Chem.Mol | None) -> str | None:
    if mol is None:
        return None
    key = Chem.MolToInchiKey(mol)
    return key.split("-")[0] if key else None


def raw_inchikey_block(inchi: str, smiles: str) -> str | None:
    if inchi:
        key = Chem.InchiToInchiKey(inchi)
        if key:
            return key.split("-")[0]
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            key = Chem.MolToInchiKey(mol)
            if key:
                return key.split("-")[0]
    return None


def load_anpdb_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open(ANPDB_PATH, encoding="utf-8") as fh:
        fh.readline()
        for raw in fh:
            fields = _sep_re.split(raw.rstrip("\n"))
            inchi_idx = next(
                (i for i, v in enumerate(fields) if v.startswith("InChI=")), -1
            )
            inchi = fields[inchi_idx] if inchi_idx >= 0 else ""
            smiles = ""
            if inchi_idx > 0:
                cand = fields[inchi_idx - 1]
                if not cand.startswith("InChI="):
                    smiles = cand
            rows.append(
                {
                    "molecule_id": fields[0] if fields else "",
                    "mol_name": fields[1] if len(fields) > 1 else "",
                    "smiles": smiles,
                    "inchi": inchi,
                }
            )
    return rows


def iter_chembl_blocks() -> Iterator[str]:
    with gzip.open(CHEMBL_PATH, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("standard_inchi_key")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > idx:
                key = parts[idx]
                if key:
                    yield key.split("-")[0]


def iter_coconut_blocks() -> Iterator[str]:
    with zipfile.ZipFile(COCONUT_ZIP) as zf:
        name = next(n for n in zf.namelist() if n.lower().endswith(".csv"))
        with zf.open(name) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(text)
            for row in reader:
                key = (row.get("standard_inchi_key") or "").strip()
                if key:
                    yield key.split("-")[0]


def cached_set(path: Path, builder) -> set[str]:
    if path.exists():
        with open(path, "rb") as fh:
            return pickle.load(fh)
    print(f"  building {path.name} (one-off)", file=sys.stderr)
    s = set(builder())
    with open(path, "wb") as fh:
        pickle.dump(s, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return s


def load_coconut_fps() -> tuple[list, list[str], list[str]]:
    if not COCONUT_FPS_CACHE.exists():
        raise SystemExit(
            f"Missing {COCONUT_FPS_CACHE}. Run scripts/run_tanimoto.py first to build it."
        )
    with open(COCONUT_FPS_CACHE, "rb") as fh:
        fps, ids, smis = pickle.load(fh)
    return fps, ids, smis


def main() -> int:
    print("Loading reference InChIKey-block sets ...", file=sys.stderr)
    chembl_blocks = cached_set(CHEMBL_BLOCKS_CACHE, iter_chembl_blocks)
    coconut_blocks = cached_set(COCONUT_BLOCKS_CACHE, iter_coconut_blocks)
    print(f"  ChEMBL blocks: {len(chembl_blocks):,}", file=sys.stderr)
    print(f"  COCONUT blocks: {len(coconut_blocks):,}", file=sys.stderr)

    print(f"Loading ANPDB", file=sys.stderr)
    rows = load_anpdb_rows()
    print(f"  {len(rows):,} rows", file=sys.stderr)

    print("Standardizing ANPDB structures ...", file=sys.stderr)
    raw_blocks: list[str | None] = []
    std_blocks: list[str | None] = []
    std_smiles: list[str] = []
    std_mols: list[Chem.Mol | None] = []
    n_std_failed = 0
    for i, r in enumerate(rows):
        raw_blocks.append(raw_inchikey_block(r["inchi"], r["smiles"]))
        mol = standardize(r["smiles"], r["inchi"])
        std_mols.append(mol)
        if mol is None:
            n_std_failed += 1
            std_blocks.append(None)
            std_smiles.append("")
        else:
            std_blocks.append(inchikey_block(mol))
            std_smiles.append(Chem.MolToSmiles(mol))
        if (i + 1) % 2500 == 0:
            print(f"  standardized {i + 1:,}/{len(rows):,}", file=sys.stderr)
    print(f"  standardization failures: {n_std_failed:,}", file=sys.stderr)

    df = pd.DataFrame(rows)
    df["inchikey_block_raw"] = raw_blocks
    df["inchikey_block_std"] = std_blocks
    df["std_smiles"] = std_smiles

    df["in_chembl_raw"] = df["inchikey_block_raw"].map(
        lambda b: bool(b) and b in chembl_blocks
    )
    df["in_coconut_raw"] = df["inchikey_block_raw"].map(
        lambda b: bool(b) and b in coconut_blocks
    )
    df["novel_raw"] = (
        df["inchikey_block_raw"].notna()
        & ~df["in_chembl_raw"]
        & ~df["in_coconut_raw"]
    )
    df["in_chembl_std"] = df["inchikey_block_std"].map(
        lambda b: bool(b) and b in chembl_blocks
    )
    df["in_coconut_std"] = df["inchikey_block_std"].map(
        lambda b: bool(b) and b in coconut_blocks
    )
    df["novel_std"] = (
        df["inchikey_block_std"].notna()
        & ~df["in_chembl_std"]
        & ~df["in_coconut_std"]
    )

    n_total = len(df)
    n_raw_novel = int(df["novel_raw"].sum())
    n_std_novel = int(df["novel_std"].sum())
    n_rescued = int((df["novel_raw"] & ~df["novel_std"]).sum())
    n_skel_changed = int(
        (df["inchikey_block_raw"].fillna("") != df["inchikey_block_std"].fillna("")).sum()
    )
    print(
        f"\nPass-1:"
        f"\n  raw novel               : {n_raw_novel:,}"
        f"\n  standardized novel      : {n_std_novel:,}"
        f"\n  rescued by salt-strip   : {n_rescued:,}"
        f"\n  skeletons that changed  : {n_skel_changed:,}",
        file=sys.stderr,
    )

    # Pass 2: Tanimoto on the standardized novel residual
    print("\nLoading cached COCONUT fingerprints ...", file=sys.stderr)
    coconut_fps, coconut_ids, coconut_smis = load_coconut_fps()
    print(f"  {len(coconut_fps):,} FPs loaded", file=sys.stderr)

    print("Computing standardized parent fingerprints for novel residual ...", file=sys.stderr)
    novel_idx = df.index[df["novel_std"]].tolist()
    q_fps = []
    q_keep = []
    for idx in novel_idx:
        mol = std_mols[idx]
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_BITS)
        q_fps.append(fp)
        q_keep.append(idx)
    print(f"  {len(q_fps):,} query FPs", file=sys.stderr)

    df["max_tanimoto"] = 0.0
    df["nearest_coconut_id"] = ""
    df["nearest_coconut_smiles"] = ""
    t0 = time.time()
    for k, (idx, qfp) in enumerate(zip(q_keep, q_fps)):
        sims = DataStructs.BulkTanimotoSimilarity(qfp, coconut_fps)
        best_i = 0
        best_s = sims[0]
        for j in range(1, len(sims)):
            if sims[j] > best_s:
                best_s = sims[j]
                best_i = j
        df.at[idx, "max_tanimoto"] = best_s
        df.at[idx, "nearest_coconut_id"] = coconut_ids[best_i]
        df.at[idx, "nearest_coconut_smiles"] = coconut_smis[best_i]
        if (k + 1) % 100 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (k + 1) * (len(q_fps) - k - 1)
            print(f"  Tanimoto: {k + 1}/{len(q_fps)}  elapsed={elapsed:,.0f}s  eta={eta:,.0f}s",
                  file=sys.stderr)

    df["has_close_analogue_std"] = df["novel_std"] & (df["max_tanimoto"] >= THRESHOLD)
    df["truly_novel_std"] = df["novel_std"] & (df["max_tanimoto"] < THRESHOLD)

    n_std_close = int(df["has_close_analogue_std"].sum())
    n_std_truly = int(df["truly_novel_std"].sum())
    print(
        f"\nPass-2 on standardized residual (threshold {THRESHOLD}):"
        f"\n  close analogue in COCONUT: {n_std_close:,}"
        f"\n  truly novel              : {n_std_truly:,}",
        file=sys.stderr,
    )

    # Outputs
    annotated = OUT / "anpdb_annotated_std.csv"
    novel_p1 = OUT / "anpdb_novel_std.csv"
    novel_p2 = OUT / "anpdb_truly_novel_std.csv"
    ids_p1 = OUT / "novel_molecule_ids_std.txt"
    ids_p2 = OUT / "truly_novel_molecule_ids_std.txt"
    summary = OUT / "summary_std.txt"

    df.to_csv(annotated, index=False)
    df_p1 = df[df["novel_std"]]
    df_p2 = df[df["truly_novel_std"]]
    df_p1.to_csv(novel_p1, index=False)
    df_p2.to_csv(novel_p2, index=False)
    with open(ids_p1, "w") as fh:
        for mid in df_p1["molecule_id"]:
            fh.write(f"{mid}\n")
    with open(ids_p2, "w") as fh:
        for mid in df_p2["molecule_id"]:
            fh.write(f"{mid}\n")

    # Tanimoto distribution on the standardized residual
    bins = [0.0, 0.50, 0.70, 0.80, 0.85, 0.90, 0.95, 1.001]
    labels = ["<0.50", "0.50-0.70", "0.70-0.80", "0.80-0.85",
              "0.85-0.90", "0.90-0.95", "0.95-1.00"]
    sub = df[df["novel_std"]]
    counts = pd.cut(sub["max_tanimoto"], bins=bins, right=False, labels=labels).value_counts().sort_index()

    lines = [
        f"ANPDB rows                       : {n_total:,}",
        f"Standardization failures         : {n_std_failed:,}",
        f"Rows where InChIKey block changed: {n_skel_changed:,}",
        "",
        f"Pass-1 raw novel                 : {n_raw_novel:,}",
        f"Pass-1 standardized novel        : {n_std_novel:,}",
        f"  rescued by salt-strip          : {n_rescued:,}",
        "",
        f"Pass-2 (Tanimoto >= {THRESHOLD}) on std residual:",
        f"  close analogue in COCONUT      : {n_std_close:,}",
        f"  truly novel                    : {n_std_truly:,}",
        "",
        "Distribution of max Tanimoto on std novel:",
    ]
    for lbl, c in counts.items():
        lines.append(f"  {lbl:<11}{c:>5}")
    lines += [
        "",
        f"Outputs:",
        f"  {annotated}",
        f"  {novel_p1}",
        f"  {novel_p2}",
        f"  {ids_p1}",
        f"  {ids_p2}  (TIGHTENED novel set)",
    ]
    summary.write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
