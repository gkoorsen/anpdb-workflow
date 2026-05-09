"""Stage 1: Generate Morgan ECFP4 fingerprint arrays for UMAP chemical space.

Produces numpy .npy arrays for:
  - 1,012 novel ANPDB compounds
  - Random 5,000-compound COCONUT sample (for background)

Output: output/chemspace/novel_fps.npy, coconut_sample_fps.npy, labels.csv
"""

from __future__ import annotations

import csv
import random
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
csv.field_size_limit(1 << 30)

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"
COCONUT_ZIP = ROOT / "data" / "coconut_csv.zip"
OUT = ROOT / "output" / "chemspace"
OUT.mkdir(parents=True, exist_ok=True)

MORGAN_RADIUS = 2
MORGAN_BITS = 2048
COCONUT_SAMPLE_SIZE = 5000
SEED = 42


def fp_to_array(fp):
    arr = np.zeros(MORGAN_BITS, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def compute_fps(smiles_list):
    fps = []
    valid_idx = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_BITS)
        fps.append(fp_to_array(fp))
        valid_idx.append(i)
    return np.array(fps, dtype=np.uint8), valid_idx


def main() -> int:
    random.seed(SEED)
    np.random.seed(SEED)

    print("Loading novel compounds ...")
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    novel_smiles = [r["std_smiles"] or r["smiles"] for _, r in novel.iterrows()]
    novel_ids = novel["molecule_id"].tolist()
    print(f"  {len(novel_smiles):,} novel compounds")

    print("Computing novel fingerprints ...")
    novel_fps, novel_valid = compute_fps(novel_smiles)
    print(f"  {len(novel_fps):,} valid FPs")

    print(f"\nLoading COCONUT sample ({COCONUT_SAMPLE_SIZE:,} compounds) ...")
    coconut_smiles = []
    with zipfile.ZipFile(COCONUT_ZIP) as zf:
        inner = [n for n in zf.namelist() if n.endswith(".csv")][0]
        with zf.open(inner) as raw:
            import io
            text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            reader = csv.DictReader(text)
            all_smiles = []
            for row in reader:
                smi = row.get("canonical_smiles") or row.get("smiles") or ""
                if smi:
                    all_smiles.append(smi)
            print(f"  total COCONUT SMILES: {len(all_smiles):,}")
            coconut_smiles = random.sample(all_smiles, min(COCONUT_SAMPLE_SIZE, len(all_smiles)))

    print(f"Computing COCONUT fingerprints ({len(coconut_smiles):,}) ...")
    coconut_fps, coconut_valid = compute_fps(coconut_smiles)
    print(f"  {len(coconut_fps):,} valid FPs")

    np.save(OUT / "novel_fps.npy", novel_fps)
    np.save(OUT / "coconut_sample_fps.npy", coconut_fps)

    labels = (
        [{"idx": i, "source": "Novel ANPDB", "compound_id": novel_ids[novel_valid[i]]}
         for i in range(len(novel_fps))]
        + [{"idx": len(novel_fps) + i, "source": "COCONUT", "compound_id": ""}
           for i in range(len(coconut_fps))]
    )
    pd.DataFrame(labels).to_csv(OUT / "labels.csv", index=False)

    print(f"\nSaved to {OUT}/")
    print(f"  novel_fps.npy: {novel_fps.shape}")
    print(f"  coconut_sample_fps.npy: {coconut_fps.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
