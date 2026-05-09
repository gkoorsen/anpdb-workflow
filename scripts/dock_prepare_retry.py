"""Retry the 5 failed ligand preparations after stripping residual Na/K salts."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
LIGANDS = ROOT / "output" / "docking" / "ligands"
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"

PY_AUTODOCK = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/python"
MK_PREP = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/mk_prepare_ligand.py"

FAILED = ["Mol_03896", "Mol_03898", "Mol_16142", "Mol_16144", "Mol_16146"]


def strip_metals(smi: str) -> str:
    # Replace coordinated Na/K to give the deprotonated parent
    smi = re.sub(r"\[O\]\[Na\]", "[O-]", smi)
    smi = re.sub(r"\[O\]\[K\]",  "[O-]", smi)
    smi = re.sub(r"\[Na\]\[O\]", "[O-]", smi)
    smi = re.sub(r"\[K\]\[O\]",  "[O-]", smi)
    return smi


def ligand_to_pdbqt(smi: str, cid: str, out_pdbqt: Path) -> bool:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    # Neutralise charges to give a vacuum-friendly molecule for docking
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) < 0:
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) < 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    sdf_path = out_pdbqt.with_suffix(".sdf")
    Chem.SDWriter(str(sdf_path)).write(mol)
    cmd = [PY_AUTODOCK, MK_PREP, "-i", str(sdf_path), "-o", str(out_pdbqt)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    sdf_path.unlink(missing_ok=True)
    return out_pdbqt.exists() and out_pdbqt.stat().st_size > 0


def main() -> int:
    n = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    smi_map = dict(zip(n["molecule_id"], n["std_smiles"]))
    for cid in FAILED:
        smi_in = smi_map.get(cid, "")
        if not smi_in:
            print(f"{cid}: no SMILES", file=sys.stderr)
            continue
        smi_clean = strip_metals(smi_in)
        out = LIGANDS / f"{cid}.pdbqt"
        ok = ligand_to_pdbqt(smi_clean, cid, out)
        print(f"{cid}: {'OK' if ok else 'STILL FAILED'}  ({out.stat().st_size if out.exists() else 0} bytes)  smiles={smi_clean[:60]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
