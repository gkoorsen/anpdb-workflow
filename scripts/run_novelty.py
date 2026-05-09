"""Flag ANPDB compounds whose structures are absent from COCONUT and ChEMBL.

Match strategy: exact match on the InChIKey first block (the 14-char skeleton
hash). This treats stereoisomers / tautomers as the same compound, which is the
right default for "is this scaffold/connectivity already known?".

Inputs (under data/):
  - ANPDB.csv                 tab-delimited, has `smiles` and `inchi`
  - coconut_csv.zip           contains a CSV with canonical_smiles + identifier
  - chembl_36_chemreps.txt.gz tab-delimited: chembl_id, canonical_smiles, standard_inchi, standard_inchi_key

Outputs (under output/):
  - anpdb_annotated.csv  every ANPDB row + in_coconut, in_chembl, novel flags
  - anpdb_novel.csv      subset where novel == True
  - summary.txt          counts
"""

from __future__ import annotations

import csv
import gzip
import io
import sys
import zipfile
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

csv.field_size_limit(sys.maxsize)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)


def inchikey_block(inchi: str | None, smiles: str | None) -> str | None:
    """Return the 14-char skeleton block of the InChIKey, or None if unparseable.

    Prefer the supplied InChI; fall back to SMILES.
    """
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


def load_anpdb(path: Path) -> pd.DataFrame:
    """Parse ANPDB.csv — 4-space-separated, with empty cells *elided* (not blanked).

    Because empty cells are dropped, column positions are not stable across rows.
    Locate `inchi` by its 'InChI=' prefix; treat the field immediately before it
    as `smiles` if that field doesn't itself start with 'InChI='. molecule_id and
    mol_name are stable at indices 0 and 1.
    """
    import re

    sep_re = re.compile(r" {4,}")
    rows: list[dict[str, str]] = []
    with open(path, encoding="utf-8") as fh:
        fh.readline()  # skip header — positions in source are unreliable
        for raw in fh:
            fields = sep_re.split(raw.rstrip("\n"))
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
                    "raw_fields": "␟".join(fields),  # keep originals for traceability
                }
            )
    df = pd.DataFrame(rows)
    df["inchikey_block"] = [
        inchikey_block(i or None, s or None)
        for i, s in zip(df["inchi"], df["smiles"])
    ]
    return df


def iter_chembl_blocks(path: Path) -> Iterator[str]:
    """Yield InChIKey first-blocks from chembl_chemreps.txt.gz."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index("standard_inchi_key")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > idx:
                key = parts[idx]
                if key:
                    yield key.split("-")[0]


def iter_coconut_blocks(zip_path: Path) -> Iterator[str]:
    """Yield InChIKey first-blocks from the COCONUT CSV inside the zip.

    Tries an `inchikey` column first; falls back to deriving from SMILES.
    """
    with zipfile.ZipFile(zip_path) as zf:
        csv_names = [n for n in zf.namelist() if n.lower().endswith(".csv")]
        if not csv_names:
            raise RuntimeError(f"No CSV inside {zip_path}")
        name = csv_names[0]
        print(f"  reading {name} from zip ...", file=sys.stderr)
        with zf.open(name) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
            reader = csv.DictReader(text)
            cols = {c.lower(): c for c in (reader.fieldnames or [])}
            ik_col = next((cols[c] for c in cols if "inchikey" in c or "inchi_key" in c), None)
            smi_col = next((cols[c] for c in cols if c == "canonical_smiles" or c == "smiles"), None)
            if ik_col is None and smi_col is None:
                raise RuntimeError(f"COCONUT CSV has no InChIKey or SMILES column. cols={list(cols)}")
            print(f"  COCONUT inchikey col={ik_col!r} smiles fallback={smi_col!r}", file=sys.stderr)
            for i, row in enumerate(reader):
                if i and i % 100_000 == 0:
                    print(f"  COCONUT rows read: {i:,}", file=sys.stderr)
                key = (row.get(ik_col) or "").strip() if ik_col else ""
                if not key and smi_col:
                    smi = (row.get(smi_col) or "").strip()
                    if smi:
                        mol = Chem.MolFromSmiles(smi)
                        if mol is not None:
                            key = Chem.MolToInchiKey(mol) or ""
                if key:
                    yield key.split("-")[0]


def main() -> int:
    anpdb_path = DATA / "ANPDB.csv"
    coconut_path = DATA / "coconut_csv.zip"
    chembl_path = DATA / "chembl_36_chemreps.txt.gz"

    print(f"Loading ANPDB from {anpdb_path}", file=sys.stderr)
    anpdb = load_anpdb(anpdb_path)
    n_total = len(anpdb)
    n_keyed = anpdb["inchikey_block"].notna().sum()
    print(f"  ANPDB rows: {n_total:,}  with InChIKey: {n_keyed:,}", file=sys.stderr)

    print("Loading ChEMBL InChIKey blocks ...", file=sys.stderr)
    chembl_blocks: set[str] = set(iter_chembl_blocks(chembl_path))
    print(f"  ChEMBL distinct blocks: {len(chembl_blocks):,}", file=sys.stderr)

    print("Loading COCONUT InChIKey blocks ...", file=sys.stderr)
    coconut_blocks: set[str] = set(iter_coconut_blocks(coconut_path))
    print(f"  COCONUT distinct blocks: {len(coconut_blocks):,}", file=sys.stderr)

    anpdb["in_chembl"] = anpdb["inchikey_block"].map(
        lambda b: bool(b) and b in chembl_blocks
    )
    anpdb["in_coconut"] = anpdb["inchikey_block"].map(
        lambda b: bool(b) and b in coconut_blocks
    )
    anpdb["novel"] = (
        anpdb["inchikey_block"].notna()
        & ~anpdb["in_chembl"]
        & ~anpdb["in_coconut"]
    )

    annotated = OUT / "anpdb_annotated.csv"
    novel = OUT / "anpdb_novel.csv"
    summary = OUT / "summary.txt"

    anpdb.to_csv(annotated, index=False)
    anpdb[anpdb["novel"]].to_csv(novel, index=False)

    n_novel = int(anpdb["novel"].sum())
    n_in_coc = int(anpdb["in_coconut"].sum())
    n_in_chembl = int(anpdb["in_chembl"].sum())
    n_in_either = int((anpdb["in_coconut"] | anpdb["in_chembl"]).sum())
    n_in_both = int((anpdb["in_coconut"] & anpdb["in_chembl"]).sum())
    n_unkeyed = int(anpdb["inchikey_block"].isna().sum())

    lines = [
        f"ANPDB rows total                    {n_total:,}",
        f"  with computable InChIKey block    {n_keyed:,}",
        f"  unparseable (excluded)            {n_unkeyed:,}",
        "",
        f"ChEMBL distinct InChIKey blocks     {len(chembl_blocks):,}",
        f"COCONUT distinct InChIKey blocks    {len(coconut_blocks):,}",
        "",
        f"ANPDB compounds in ChEMBL           {n_in_chembl:,}",
        f"ANPDB compounds in COCONUT          {n_in_coc:,}",
        f"ANPDB compounds in either           {n_in_either:,}",
        f"ANPDB compounds in both             {n_in_both:,}",
        "",
        f"ANPDB compounds NOVEL (in neither)  {n_novel:,}  "
        f"({100 * n_novel / n_total:.1f}% of total, "
        f"{100 * n_novel / max(n_keyed, 1):.1f}% of keyed)",
        "",
        f"Outputs:",
        f"  {annotated}",
        f"  {novel}",
    ]
    summary.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
