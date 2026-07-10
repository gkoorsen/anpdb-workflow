"""Diagnose ligand coordinate grafting used by the OpenMM/OpenFF MD workflow."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]

CASES = [
    ("Mol_13144", "data/md_inputs/poses/Mol_13144_7VSI_opm_oriented_out.pdbqt"),
    ("Mol_13733", "data/md_inputs/poses/Mol_13733_7VSI_opm_oriented_out.pdbqt"),
    ("Mol_15088", "data/md_inputs/poses/Mol_15088_7VSI_opm_oriented_out.pdbqt"),
    ("Mol_14056", "data/md_inputs/poses/Mol_14056_2V5Z_out.pdbqt"),
    ("Mol_16614", "data/md_inputs/poses/Mol_16614_4DJH_opm_oriented_out.pdbqt"),
]


def extract_pdbqt_model(pdbqt_path: Path, mode: int, out_pdbqt: Path) -> None:
    keep: list[str] = []
    in_mode = False
    saw = 0
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL"):
            saw += 1
            in_mode = saw == mode
            if in_mode:
                keep.append(line)
            continue
        if line.startswith("ENDMDL"):
            if in_mode:
                keep.append(line)
                break
            in_mode = False
            continue
        if in_mode:
            keep.append(line)
    if not keep:
        raise RuntimeError(f"Could not extract mode {mode} from {pdbqt_path}")
    out_pdbqt.write_text("\n".join(keep) + "\n", encoding="utf-8")


def convert_pdbqt_to_pdb(pdbqt_path: Path, out_pdb: Path, obabel: str) -> None:
    result = subprocess.run(
        [obabel, str(pdbqt_path), "-O", str(out_pdb)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout)


def heavy_pose_from_pdb(pose_pdb: Path):
    raw = Chem.MolFromPDBFile(str(pose_pdb), removeHs=False, sanitize=False)
    if raw is None:
        raise RuntimeError(f"RDKit could not read {pose_pdb}")
    editable = Chem.EditableMol(raw)
    hydrogens = sorted([atom.GetIdx() for atom in raw.GetAtoms() if atom.GetSymbol() == "H"], reverse=True)
    for atom_index in hydrogens:
        editable.RemoveAtom(atom_index)
    raw_no_h = editable.GetMol()
    Chem.SanitizeMol(raw_no_h)
    return raw_no_h


def bond_length_stats(canonical_no_h, pose_with_bonds, canonical_to_pose: dict[int, int]) -> dict[str, float | int]:
    pose_conf = pose_with_bonds.GetConformer()
    coords = []
    for canonical_index in range(canonical_no_h.GetNumAtoms()):
        pose_index = canonical_to_pose[canonical_index]
        point = pose_conf.GetAtomPosition(pose_index)
        coords.append((point.x, point.y, point.z))
    lengths = [
        math.dist(coords[bond.GetBeginAtomIdx()], coords[bond.GetEndAtomIdx()])
        for bond in canonical_no_h.GetBonds()
    ]
    bad = [length for length in lengths if length < 0.90 or length > 1.90]
    return {
        "min_bond_A": min(lengths),
        "max_bond_A": max(lengths),
        "bad_bond_count": len(bad),
    }


def mapping_rows(compound_id: str, pdbqt_path: Path, smiles: str, obabel: str, work_dir: Path) -> list[dict[str, object]]:
    mode_pdbqt = work_dir / f"{compound_id}_mode1.pdbqt"
    pose_pdb = work_dir / f"{compound_id}_mode1.pdb"
    extract_pdbqt_model(pdbqt_path, 1, mode_pdbqt)
    convert_pdbqt_to_pdb(mode_pdbqt, pose_pdb, obabel)

    canonical = Chem.MolFromSmiles(smiles)
    if canonical is None:
        raise RuntimeError(f"RDKit could not parse SMILES for {compound_id}")
    canonical = Chem.AddHs(canonical)
    AllChem.EmbedMolecule(canonical, randomSeed=42, useRandomCoords=True)
    canonical_no_h = Chem.RemoveHs(canonical, sanitize=True)
    raw_no_h = heavy_pose_from_pdb(pose_pdb)
    pose_with_bonds = AllChem.AssignBondOrdersFromTemplate(canonical_no_h, raw_no_h)

    rows: list[dict[str, object]] = []
    first_match = canonical_no_h.GetSubstructMatch(pose_with_bonds)
    pose_match = pose_with_bonds.GetSubstructMatch(canonical_no_h)
    mappings = {
        "old_md_code": {canonical_index: first_match[canonical_index] for canonical_index in range(len(first_match))},
        "invert_old_match": {canonical_index: pose_index for pose_index, canonical_index in enumerate(first_match)},
        "fixed_pose_match": {canonical_index: pose_index for canonical_index, pose_index in enumerate(pose_match)},
    }
    for mapping_name, mapping in mappings.items():
        stats = bond_length_stats(canonical_no_h, pose_with_bonds, mapping)
        rows.append(
            {
                "compound_id": compound_id,
                "mapping": mapping_name,
                "heavy_atoms": canonical_no_h.GetNumAtoms(),
                "heavy_bonds": canonical_no_h.GetNumBonds(),
                "first_match_head": ";".join(str(item) for item in first_match[:12]),
                "pose_match_head": ";".join(str(item) for item in pose_match[:12]),
                **stats,
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--out", type=Path, default=Path("results/md_analysis/ligand_graft_mapping_diagnostics.csv"))
    parser.add_argument("--smiles-table", type=Path, default=Path("data/md_inputs/anpdb_truly_novel_std.csv"))
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    smiles_table = pd.read_csv(ROOT / args.smiles_table).set_index("molecule_id")
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory() as tmp:
        work_dir = Path(tmp)
        for compound_id, pdbqt_rel in CASES:
            rows.extend(mapping_rows(compound_id, ROOT / pdbqt_rel, smiles_table.loc[compound_id, "std_smiles"], args.obabel, work_dir))

    out_path = args.out if args.out.is_absolute() else ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
