"""Build curated Amber production input systems.

This builder is intentionally conservative. It creates the CYP1B1 heme system
using the Shahrokh/Orendt/Yost/Cheatham P450 IC6 heme parameters bundled with
AmberTools. It also creates an explicitly versioned MAO-B CYF residue for
8alpha-S-cysteinyl-FAD from a capped Cys397-FAD AM1-BCC model.

Run from the repository root with the Docking environment:

  conda run -n Docking python scripts/md_build_amber_systems.py --target cyp1b1 --force
  conda run -n Docking python scripts/md_build_amber_systems.py --target maob --force
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "md_inputs"
AMBER_OUT = DATA / "amber_systems"
COFACTOR_AMBER = DATA / "cofactors" / "amber"
BUILD_ROOT = ROOT / "runs" / "amber_system_prep"


def find_amberhome() -> Path:
    if os.environ.get("AMBERHOME"):
        return Path(os.environ["AMBERHOME"])
    tleap = shutil.which("tleap")
    if tleap:
        return Path(tleap).resolve().parents[1]
    return Path(sys.prefix)


AMBERHOME = find_amberhome()
SHAHROKH_IC6 = AMBERHOME / "dat" / "contrib" / "Shahrokh_heme" / "IC6"


@dataclass(frozen=True)
class LigandSpec:
    compound_id: str
    pdb_id: str
    pose_pdbqt: Path
    net_charge: int


def fail(message: str) -> None:
    raise SystemExit(message)


def run(cmd: list[str], cwd: Path, log_name: str) -> subprocess.CompletedProcess[str]:
    cwd.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False)
    (cwd / f"{log_name}.stdout").write_text(result.stdout)
    (cwd / f"{log_name}.stderr").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}) in {cwd}: {' '.join(cmd)}\n"
            f"stderr: {result.stderr[-3000:]}"
        )
    return result


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        fail(f"Required tool not found on PATH: {name}. Run through `conda run -n Docking ...`.")
    return path


def require_file(path: Path) -> Path:
    if not path.exists():
        fail(f"Missing required input: {path}")
    return path


def require_shahrokh_ic6() -> None:
    for name in ("HEM.mol2", "CYP.mol2", "IC6.frcmod"):
        require_file(SHAHROKH_IC6 / name)


def smiles_for(compound_id: str) -> str:
    table = require_file(DATA / "anpdb_truly_novel_std.csv")
    with table.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("molecule_id") == compound_id:
                smiles = row.get("std_smiles", "").strip()
                if not smiles:
                    fail(f"No std_smiles for {compound_id} in {table}")
                return smiles
    fail(f"No row for {compound_id} in {table}")


def extract_pdbqt_model(pdbqt_path: Path, out_pdbqt: Path, mode: int = 1) -> None:
    keep: list[str] = []
    current = 0
    in_model = False
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL"):
            current += 1
            in_model = current == mode
            if in_model:
                keep.append(line)
            continue
        if line.startswith("ENDMDL"):
            if in_model:
                keep.append(line)
                break
            in_model = False
            continue
        if in_model:
            keep.append(line)
    if not keep:
        fail(f"Could not extract MODEL {mode} from {pdbqt_path}")
    out_pdbqt.write_text("\n".join(keep) + "\n")


def pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path, work_dir: Path) -> None:
    run([require_tool("obabel"), str(pdbqt_path), "-O", str(pdb_path)], work_dir, "obabel_pdbqt_to_pdb")
    require_file(pdb_path)


def graft_pose_to_sdf(pose_pdb: Path, smiles: str, out_sdf: Path) -> None:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    canonical = Chem.MolFromSmiles(smiles)
    if canonical is None:
        fail(f"RDKit failed to parse SMILES: {smiles}")
    canonical = Chem.AddHs(canonical)
    if AllChem.EmbedMolecule(canonical, randomSeed=42) < 0:
        AllChem.EmbedMolecule(canonical, randomSeed=42, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(canonical, maxIters=200)

    raw = Chem.MolFromPDBFile(str(pose_pdb), removeHs=False, sanitize=False)
    if raw is None:
        fail(f"RDKit failed to read pose PDB: {pose_pdb}")
    editable = Chem.EditableMol(raw)
    for idx in sorted((a.GetIdx() for a in raw.GetAtoms() if a.GetSymbol() == "H"), reverse=True):
        editable.RemoveAtom(idx)
    raw_no_h = editable.GetMol()
    Chem.SanitizeMol(raw_no_h)
    canonical_no_h = Chem.RemoveHs(canonical, sanitize=True)

    if canonical_no_h.GetNumAtoms() != raw_no_h.GetNumAtoms():
        fail(
            "Pose/canonical heavy-atom count mismatch: "
            f"{raw_no_h.GetNumAtoms()} vs {canonical_no_h.GetNumAtoms()}"
        )

    posed_with_bonds = AllChem.AssignBondOrdersFromTemplate(canonical_no_h, raw_no_h)
    match = posed_with_bonds.GetSubstructMatch(canonical_no_h)
    if not match or len(match) != canonical_no_h.GetNumAtoms():
        reverse = canonical_no_h.GetSubstructMatch(posed_with_bonds)
        if reverse and len(reverse) == canonical_no_h.GetNumAtoms():
            inverse = [0] * len(reverse)
            for pose_idx, canonical_idx in enumerate(reverse):
                inverse[canonical_idx] = pose_idx
            match = tuple(inverse)
    if not match or len(match) != canonical_no_h.GetNumAtoms():
        fail("Could not map docked pose atoms onto canonical ligand graph")

    pose_conf = posed_with_bonds.GetConformer()
    canonical_conf = canonical.GetConformer()
    for canonical_idx, pose_idx in enumerate(match):
        canonical_conf.SetAtomPosition(canonical_idx, pose_conf.GetAtomPosition(pose_idx))

    canonical_heavy = Chem.RemoveHs(canonical)
    canonical_with_h = Chem.AddHs(canonical_heavy, addCoords=True)
    Chem.SanitizeMol(canonical_with_h)
    conf = canonical_with_h.GetConformer()
    long_bonds: list[str] = []
    for bond in canonical_with_h.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() == 1 or b.GetAtomicNum() == 1:
            continue
        pa = conf.GetAtomPosition(a.GetIdx())
        pb = conf.GetAtomPosition(b.GetIdx())
        distance = math.dist((pa.x, pa.y, pa.z), (pb.x, pb.y, pb.z))
        if distance > 2.1:
            long_bonds.append(f"{a.GetSymbol()}{a.GetIdx()}-{b.GetSymbol()}{b.GetIdx()}={distance:.2f}A")
    if long_bonds:
        fail("Mapped ligand pose has implausible heavy-atom bond lengths: " + ", ".join(long_bonds[:12]))

    writer = Chem.SDWriter(str(out_sdf))
    writer.write(canonical_with_h)
    writer.close()
    require_file(out_sdf)


def parameterize_ligand(spec: LigandSpec, work_dir: Path, force: bool) -> tuple[Path, Path]:
    mode_pdbqt = work_dir / f"{spec.compound_id}_mode1.pdbqt"
    pose_pdb = work_dir / f"{spec.compound_id}_pose.pdb"
    ligand_sdf = work_dir / f"{spec.compound_id}.sdf"
    ligand_mol2 = work_dir / f"{spec.compound_id}_gaff2.mol2"
    ligand_frcmod = work_dir / f"{spec.compound_id}_gaff2.frcmod"

    if not force and ligand_mol2.exists() and ligand_frcmod.exists():
        return ligand_mol2, ligand_frcmod

    extract_pdbqt_model(require_file(spec.pose_pdbqt), mode_pdbqt)
    pdbqt_to_pdb(mode_pdbqt, pose_pdb, work_dir)
    graft_pose_to_sdf(pose_pdb, smiles_for(spec.compound_id), ligand_sdf)

    run(
        [
            require_tool("antechamber"),
            "-i",
            str(ligand_sdf),
            "-fi",
            "sdf",
            "-o",
            str(ligand_mol2),
            "-fo",
            "mol2",
            "-rn",
            "LIG",
            "-at",
            "gaff2",
            "-c",
            "bcc",
            "-nc",
            str(spec.net_charge),
            "-s",
            "2",
            "-pf",
            "y",
        ],
        work_dir,
        "antechamber_ligand",
    )
    run(
        [
            require_tool("parmchk2"),
            "-i",
            str(ligand_mol2),
            "-f",
            "mol2",
            "-o",
            str(ligand_frcmod),
            "-s",
            "gaff2",
        ],
        work_dir,
        "parmchk2_ligand",
    )
    return require_file(ligand_mol2), require_file(ligand_frcmod)


def prepare_cyp1b1_receptor(work_dir: Path) -> Path:
    source = require_file(DATA / "receptors" / "4I8V_chainA_heme_prepared.pdb")
    out = work_dir / "cyp1b1_cyp457_receptor.pdb"
    lines: list[str] = []
    previous_c: tuple[float, float, float] | None = None
    cyp_n: tuple[float, float, float] | None = None
    cyp_ca: tuple[float, float, float] | None = None
    cyp_n_insert_at: int | None = None
    for line in source.read_text().splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM", "TER"}:
            continue
        chain = line[21:22]
        if chain != "A":
            continue
        if record in {"ATOM", "HETATM"}:
            resname = line[17:20]
            resseq = int(line[22:26])
            atom_name = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if resname == "LYS" and resseq == 456 and atom_name == "C":
                previous_c = xyz
            if resname == "CYS" and resseq == 457:
                line = f"{line[:17]}CYP{line[20:]}"
                if atom_name == "N":
                    cyp_n = xyz
                    cyp_n_insert_at = len(lines) + 1
                elif atom_name == "CA":
                    cyp_ca = xyz
        lines.append(line)
    if previous_c is None or cyp_n is None or cyp_ca is None or cyp_n_insert_at is None:
        fail("Could not locate LYS456 C and CYP457 N/CA needed to place CYP457 backbone H")

    def unit_vector(vector: tuple[float, float, float]) -> tuple[float, float, float]:
        length = math.sqrt(sum(item * item for item in vector))
        if length == 0:
            fail("Zero-length vector while placing CYP457 backbone H")
        return tuple(item / length for item in vector)

    away_from_ca = unit_vector(tuple(cyp_n[i] - cyp_ca[i] for i in range(3)))
    away_from_prev_c = unit_vector(tuple(cyp_n[i] - previous_c[i] for i in range(3)))
    direction = unit_vector(tuple(away_from_ca[i] + away_from_prev_c[i] for i in range(3)))
    h_xyz = tuple(cyp_n[i] + 1.01 * direction[i] for i in range(3))
    h_line = (
        f"ATOM  {99999:5d}  H   CYP A{457:4d}    "
        f"{h_xyz[0]:8.3f}{h_xyz[1]:8.3f}{h_xyz[2]:8.3f}"
        "  1.00  0.00           H  "
    )
    lines.insert(cyp_n_insert_at, h_line)
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    return out


def parse_box_volume(inpcrd: Path) -> float:
    last = ""
    for line in inpcrd.read_text().splitlines():
        if line.strip():
            last = line
    values = [float(item) for item in last.split()]
    if len(values) < 3:
        fail(f"Could not parse periodic box from {inpcrd}")
    return values[0] * values[1] * values[2]


def salt_pairs_for_volume(volume_a3: float, molar: float) -> int:
    return max(0, int(round(molar * 0.000602214076 * volume_a3)))


def renumber_pdb_lines(lines: list[str]) -> list[str]:
    out: list[str] = []
    serial = 1
    for line in lines:
        if line[:6].strip() in {"ATOM", "HETATM"}:
            out.append(f"{line[:6]}{serial:5d}{line[11:]}")
            serial += 1
        else:
            out.append(line)
    return out


def write_cyp1b1_tleap(
    path: Path,
    receptor_pdb: Path,
    ligand_mol2: Path,
    ligand_frcmod: Path,
    out_prefix: Path,
    salt_pairs: int,
) -> None:
    lines = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
        f'HEM = loadmol2 "{SHAHROKH_IC6 / "HEM.mol2"}"',
        f'CYP = loadmol2 "{SHAHROKH_IC6 / "CYP.mol2"}"',
        f'loadamberparams "{SHAHROKH_IC6 / "IC6.frcmod"}"',
        f'LIG = loadmol2 "{ligand_mol2}"',
        f'loadamberparams "{ligand_frcmod}"',
        f'REC = loadpdb "{receptor_pdb}"',
        "bond REC.456.C REC.457.N",
        "bond REC.457.C REC.458.N",
        "bond REC.457.SG REC.513.FE",
        "COM = combine { REC LIG }",
        "check COM",
        "charge COM",
        "solvateBox COM TIP3PBOX 10.0",
        "addIonsRand COM Na+ 0",
        "addIonsRand COM Cl- 0",
    ]
    if salt_pairs > 0:
        lines.extend(
            [
                f"addIonsRand COM Na+ {salt_pairs}",
                f"addIonsRand COM Cl- {salt_pairs}",
            ]
        )
    lines.extend(
        [
            "check COM",
            f'saveamberparm COM "{out_prefix.with_suffix(".prmtop")}" "{out_prefix.with_suffix(".inpcrd")}"',
            f'savepdb COM "{out_prefix.with_suffix(".pdb")}"',
            "quit",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def tleap_build_cyp1b1(force: bool) -> Path:
    final_prmtop = AMBER_OUT / "cyp1b1_mol11315.prmtop"
    final_inpcrd = AMBER_OUT / "cyp1b1_mol11315.inpcrd"
    if not force and final_prmtop.exists() and final_inpcrd.exists():
        return final_prmtop

    require_shahrokh_ic6()
    AMBER_OUT.mkdir(parents=True, exist_ok=True)
    work_dir = BUILD_ROOT / "cyp1b1_mol11315"
    work_dir.mkdir(parents=True, exist_ok=True)

    ligand = LigandSpec(
        compound_id="Mol_11315",
        pdb_id="4I8V",
        pose_pdbqt=DATA / "poses" / "Mol_11315_4I8V_out.pdbqt",
        net_charge=0,
    )
    ligand_mol2, ligand_frcmod = parameterize_ligand(ligand, work_dir, force=force)
    receptor_pdb = prepare_cyp1b1_receptor(work_dir)

    estimate_prefix = work_dir / "cyp1b1_mol11315_neutral_estimate"
    estimate_tleap = work_dir / "tleap_cyp1b1_estimate.in"
    write_cyp1b1_tleap(estimate_tleap, receptor_pdb, ligand_mol2, ligand_frcmod, estimate_prefix, 0)
    run([require_tool("tleap"), "-f", str(estimate_tleap)], work_dir, "tleap_cyp1b1_estimate")
    volume_a3 = parse_box_volume(estimate_prefix.with_suffix(".inpcrd"))
    salt_pairs = salt_pairs_for_volume(volume_a3, 0.15)

    final_prefix = AMBER_OUT / "cyp1b1_mol11315"
    final_tleap = work_dir / "tleap_cyp1b1_final.in"
    write_cyp1b1_tleap(final_tleap, receptor_pdb, ligand_mol2, ligand_frcmod, final_prefix, salt_pairs)
    run([require_tool("tleap"), "-f", str(final_tleap)], work_dir, "tleap_cyp1b1_final")

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "system": "cyp1b1_mol11315",
        "target": "CYP1B1",
        "compound_id": "Mol_11315",
        "receptor": str(DATA / "receptors" / "4I8V_chainA_heme_prepared.pdb"),
        "ligand_pose": str(ligand.pose_pdbqt),
        "ligand_net_charge": ligand.net_charge,
        "protein_force_field": "Amber ff14SB",
        "ligand_force_field": "GAFF2 with AM1-BCC charges from antechamber",
        "heme_parameters": {
            "source": "AmberTools dat/contrib/Shahrokh_heme/IC6",
            "state": "Ferric penta-coordinate high-spin P450 heme, proximal cysteine ligand",
            "files": [
                str(SHAHROKH_IC6 / "HEM.mol2"),
                str(SHAHROKH_IC6 / "CYP.mol2"),
                str(SHAHROKH_IC6 / "IC6.frcmod"),
            ],
            "bond": "CYP457 SG - HEM513 FE in Leap numbering; source PDB residue HEM A 601",
        },
        "water_model": "TIP3PBOX",
        "padding_angstrom": 10.0,
        "target_salt_molar": 0.15,
        "estimated_neutral_box_volume_a3": volume_a3,
        "added_salt_pairs": salt_pairs,
        "outputs": {
            "prmtop": str(final_prmtop),
            "inpcrd": str(final_inpcrd),
        },
        "work_dir": str(work_dir),
    }
    (work_dir / "cyp1b1_mol11315_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return require_file(final_prmtop)


def write_maob_cyf_capped_heavy(work_dir: Path) -> Path:
    source = require_file(DATA / "receptors" / "2V5Z_chainA_fad_prepared.pdb")
    out = work_dir / "cyf_capped_heavy.pdb"
    records: list[tuple[int, str, str, str]] = []
    for line in source.read_text().splitlines():
        if line[:6].strip() in {"ATOM", "HETATM"} and line[21:22] == "A":
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            atom_name = line[12:16].strip()
            if resseq in {396, 397, 398} or resname == "FAD":
                records.append((resseq, resname, atom_name, line))

    keep: list[tuple[str, str, str, int]] = []
    for resseq, _resname, atom_name, line in records:
        if resseq == 396 and atom_name in {"CA", "C", "O"}:
            keep.append(({"CA": "PCA", "C": "PC", "O": "PO"}[atom_name], line, "CAP", 396))
    for resseq, resname, atom_name, line in records:
        if resseq == 397 and resname == "CYS":
            keep.append((atom_name, line, "CYF", 397))
    for _resseq, resname, atom_name, line in records:
        if resname == "FAD":
            keep.append((atom_name, line, "CYF", 397))
    for resseq, _resname, atom_name, line in records:
        if resseq == 398 and atom_name in {"N", "CA"}:
            keep.append(({"N": "NN", "CA": "NCA"}[atom_name], line, "CAP", 398))

    serial_by_name: dict[str, int] = {}
    lines: list[str] = []
    for serial, (name, line, resname, resseq) in enumerate(keep, 1):
        serial_by_name[name] = serial
        element = (line[76:78].strip() or name[0]).rjust(2)
        lines.append(
            f"HETATM{serial:5d} {name:<4s} {resname:>3s} A{resseq:4d}    "
            f"{line[30:38]}{line[38:46]}{line[46:54]}  1.00  0.00          {element}"
        )
    bonds = [
        ("PCA", "PC"),
        ("PC", "PO"),
        ("PC", "N"),
        ("N", "CA"),
        ("CA", "C"),
        ("C", "O"),
        ("CA", "CB"),
        ("CB", "SG"),
        ("SG", "C8M"),
        ("C", "NN"),
        ("NN", "NCA"),
    ]
    for atom1, atom2 in bonds:
        if atom1 in serial_by_name and atom2 in serial_by_name:
            lines.append(f"CONECT{serial_by_name[atom1]:5d}{serial_by_name[atom2]:5d}")
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    return out


def rename_cyf_cap_atoms(ac_path: Path, out_path: Path, mainchain_path: Path) -> None:
    renames = {
        1: "PCA",
        2: "PC",
        3: "PO",
        63: "NN",
        64: "NCA",
        65: "PH1",
        66: "PH2",
        67: "PH3",
        68: "HN",
        106: "NHN",
        107: "NM1",
        108: "NM2",
        109: "NM3",
    }
    lines: list[str] = []
    for line in ac_path.read_text().splitlines():
        if line.startswith("ATOM"):
            parts = line.split()
            idx = int(parts[1])
            name = renames.get(idx, parts[2])
            resname = parts[3]
            resseq = int(parts[4])
            x, y, z = map(float, parts[5:8])
            charge = float(parts[8])
            atom_type = parts[9]
            line = (
                f"ATOM {idx:6d} {name:<4s} {resname:<4s} {resseq:5d} "
                f"{x:11.3f} {y:8.3f} {z:8.3f} {charge:10.6f} {atom_type:>9s}"
            )
        lines.append(line)
    out_path.write_text("\n".join(lines) + "\n")
    mainchain_path.write_text(
        "\n".join(
            [
                "HEAD_NAME N",
                "TAIL_NAME C",
                "MAIN_CHAIN CA",
                "OMIT_NAME PCA",
                "OMIT_NAME PC",
                "OMIT_NAME PO",
                "OMIT_NAME PH1",
                "OMIT_NAME PH2",
                "OMIT_NAME PH3",
                "OMIT_NAME NN",
                "OMIT_NAME NCA",
                "OMIT_NAME NHN",
                "OMIT_NAME NM1",
                "OMIT_NAME NM2",
                "OMIT_NAME NM3",
                "PRE_HEAD_TYPE C",
                "POST_TAIL_TYPE N",
                "CHARGE -2.0",
            ]
        )
        + "\n"
    )


def write_stripped_cyf_mol2(ac_path: Path, mol2_path: Path, manifest_path: Path) -> None:
    omit = {1, 2, 3, 63, 64, 65, 66, 67, 106, 107, 108, 109}
    atoms: dict[int, dict[str, object]] = {}
    bonds: list[tuple[int, int, str]] = []
    for line in ac_path.read_text().splitlines():
        if line.startswith("ATOM"):
            parts = line.split()
            atoms[int(parts[1])] = {
                "name": parts[2],
                "type": parts[9],
                "charge": float(parts[8]),
                "xyz": tuple(float(item) for item in parts[5:8]),
            }
        elif line.startswith("BOND"):
            parts = line.split()
            bonds.append((int(parts[2]), int(parts[3]), parts[4]))

    kept = [idx for idx in sorted(atoms) if idx not in omit]
    raw_charge = sum(float(atoms[idx]["charge"]) for idx in kept)
    correction = (-2.0 - raw_charge) / len(kept)
    old_to_new = {old: new for new, old in enumerate(kept, 1)}
    protein_types = {
        "N": "N",
        "HN": "H",
        "CA": "CT",
        "C": "C",
        "O": "O",
        "CB": "CT",
        "H3": "H1",
        "H4": "H1",
        "H5": "H1",
        "SG": "S",
    }
    kept_bonds = [(a, b, order) for a, b, order in bonds if a in old_to_new and b in old_to_new]
    lines = [
        "@<TRIPOS>MOLECULE",
        "CYF",
        f"{len(kept)} {len(kept_bonds)} 1 0 0",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    for old_idx in kept:
        new_idx = old_to_new[old_idx]
        atom = atoms[old_idx]
        name = str(atom["name"])
        atom_type = protein_types.get(name, str(atom["type"]))
        x, y, z = atom["xyz"]  # type: ignore[misc]
        charge = float(atom["charge"]) + correction
        lines.append(
            f"{new_idx:7d} {name:<4s} {x:10.4f} {y:10.4f} {z:10.4f} "
            f"{atom_type:<6s} 1 CYF {charge:12.6f}"
        )
    lines.append("@<TRIPOS>BOND")
    for bond_idx, (atom1, atom2, order) in enumerate(kept_bonds, 1):
        lines.append(f"{bond_idx:6d} {old_to_new[atom1]:5d} {old_to_new[atom2]:5d} {order}")
    lines.extend(["@<TRIPOS>SUBSTRUCTURE", "     1 CYF         1 RESIDUE           0 ****  ****    0 ROOT"])
    mol2_path.write_text("\n".join(lines) + "\n")
    manifest = {
        "residue": "CYF",
        "description": "Cys397 covalently linked to oxidized FAD through SG-C8M thioether",
        "charge_model": "AmberTools antechamber GAFF2 AM1-BCC on capped Cys-FAD model",
        "target_residue_charge": -2.0,
        "raw_stripped_charge": raw_charge,
        "uniform_charge_correction_per_atom": correction,
        "stripped_cap_charge": sum(float(atoms[idx]["charge"]) for idx in omit),
        "atoms": len(kept),
        "bonds": len(kept_bonds),
        "protein_typed_atoms": sorted(protein_types),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def filter_cyf_frcmod(src: Path, out: Path) -> None:
    standard = {"N", "CT", "C", "O", "S", "H", "H1"}
    section: str | None = None
    lines: list[str] = []
    for line in src.read_text().splitlines():
        stripped = line.strip()
        if stripped in {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}:
            section = stripped
            lines.append(line)
            continue
        if not stripped or stripped.startswith("Remark"):
            lines.append(line)
            continue
        if section in {"MASS", "NONBON"}:
            if stripped.split()[0] in standard:
                continue
            lines.append(line)
            continue
        if section in {"BOND", "ANGLE", "DIHE", "IMPROPER"}:
            match = re.search(r"\s[-+]?\d+(?:\.\d+)?", line)
            param_part = line[: match.start()] if match else line
            if not re.search(r"[a-z]", param_part):
                continue
            lines.append(line)
            continue
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")


def build_maob_cyf_parameters(force: bool) -> tuple[Path, Path]:
    final_mol2 = COFACTOR_AMBER / "CYF_cys397_fad.mol2"
    final_frcmod = COFACTOR_AMBER / "CYF_cys397_fad.frcmod"
    final_manifest = COFACTOR_AMBER / "CYF_cys397_fad_manifest.json"
    if not force and final_mol2.exists() and final_frcmod.exists() and final_manifest.exists():
        return final_mol2, final_frcmod

    COFACTOR_AMBER.mkdir(parents=True, exist_ok=True)
    work_dir = BUILD_ROOT / "maob_mol14056" / "cyf_parameters"
    work_dir.mkdir(parents=True, exist_ok=True)
    heavy_pdb = write_maob_cyf_capped_heavy(work_dir)
    capped_mol2 = work_dir / "cyf_capped_h_min.mol2"
    run(
        [
            require_tool("obabel"),
            str(heavy_pdb),
            "-O",
            str(capped_mol2),
            "-h",
            "--minimize",
            "--ff",
            "GAFF",
            "--steps",
            "250",
        ],
        work_dir,
        "obabel_cyf_capped",
    )
    capped_ac = work_dir / "cyf_capped_gaff2.ac"
    run(
        [
            require_tool("antechamber"),
            "-i",
            str(capped_mol2),
            "-fi",
            "mol2",
            "-o",
            str(capped_ac),
            "-fo",
            "ac",
            "-rn",
            "CYF",
            "-at",
            "gaff2",
            "-c",
            "bcc",
            "-nc",
            "-2",
            "-s",
            "2",
            "-pf",
            "y",
        ],
        work_dir,
        "antechamber_cyf_capped",
    )
    renamed_ac = work_dir / "cyf_capped_gaff2_renamed.ac"
    rename_cyf_cap_atoms(capped_ac, renamed_ac, work_dir / "cyf.mc")
    mol2 = work_dir / "CYF.mol2"
    manifest = work_dir / "CYF_cys397_fad_manifest.json"
    write_stripped_cyf_mol2(renamed_ac, mol2, manifest)
    raw_frcmod = work_dir / "CYF.raw.frcmod"
    run([require_tool("parmchk2"), "-i", str(mol2), "-f", "mol2", "-o", str(raw_frcmod), "-s", "gaff2"], work_dir, "parmchk2_cyf")
    frcmod = work_dir / "CYF.frcmod"
    filter_cyf_frcmod(raw_frcmod, frcmod)
    shutil.copyfile(mol2, final_mol2)
    shutil.copyfile(frcmod, final_frcmod)
    shutil.copyfile(manifest, final_manifest)
    return require_file(final_mol2), require_file(final_frcmod)


def prepare_maob_receptor(work_dir: Path) -> Path:
    source = require_file(DATA / "receptors" / "2V5Z_chainA_fad_prepared.pdb")
    out = work_dir / "maob_cyf397_receptor.pdb"
    source_lines = source.read_text().splitlines()
    fad_lines = [
        f"HETATM{0:5d} {line[12:16]} CYF A{397:4d}{line[26:]}"
        for line in source_lines
        if line[:6].strip() in {"ATOM", "HETATM"} and line[21:22] == "A" and line[17:20].strip() == "FAD"
    ]
    lines: list[str] = []
    inserted = False
    for line in source_lines:
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM", "TER"}:
            continue
        if record in {"ATOM", "HETATM"}:
            if line[21:22] != "A":
                continue
            resname = line[17:20].strip()
            resseq = int(line[22:26])
            atom_name = line[12:16].strip()
            if resname == "FAD":
                continue
            if resname == "CYS" and resseq == 397:
                line = f"{line[:17]}CYF{line[20:]}"
            lines.append(line)
            if resname == "CYS" and resseq == 397 and atom_name == "SG":
                lines.extend(fad_lines)
                inserted = True
        else:
            lines.append(line)
    if not inserted:
        fail("Could not insert FAD atoms into CYF397 receptor model")
    lines = renumber_pdb_lines(lines)
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    return out


def write_maob_tleap(
    path: Path,
    receptor_pdb: Path,
    ligand_mol2: Path,
    ligand_frcmod: Path,
    cyf_mol2: Path,
    cyf_frcmod: Path,
    out_prefix: Path,
    salt_pairs: int,
) -> None:
    lines = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
        f'loadamberparams "{cyf_frcmod}"',
        f'loadamberparams "{ligand_frcmod}"',
        f'CYF = loadmol2 "{cyf_mol2}"',
        f'LIG = loadmol2 "{ligand_mol2}"',
        f'REC = loadpdb "{receptor_pdb}"',
        "bond REC.396.C REC.397.N",
        "bond REC.397.C REC.398.N",
        "COM = combine { REC LIG }",
        "check COM",
        "charge COM",
        "solvateBox COM TIP3PBOX 10.0",
        "addIonsRand COM Na+ 0",
        "addIonsRand COM Cl- 0",
    ]
    if salt_pairs > 0:
        lines.extend([f"addIonsRand COM Na+ {salt_pairs}", f"addIonsRand COM Cl- {salt_pairs}"])
    lines.extend(
        [
            "check COM",
            f'saveamberparm COM "{out_prefix.with_suffix(".prmtop")}" "{out_prefix.with_suffix(".inpcrd")}"',
            f'savepdb COM "{out_prefix.with_suffix(".pdb")}"',
            "quit",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def tleap_build_maob(force: bool) -> Path:
    final_prmtop = AMBER_OUT / "maob_mol14056.prmtop"
    final_inpcrd = AMBER_OUT / "maob_mol14056.inpcrd"
    if not force and final_prmtop.exists() and final_inpcrd.exists():
        return final_prmtop

    AMBER_OUT.mkdir(parents=True, exist_ok=True)
    work_dir = BUILD_ROOT / "maob_mol14056"
    work_dir.mkdir(parents=True, exist_ok=True)
    cyf_mol2, cyf_frcmod = build_maob_cyf_parameters(force=force)
    ligand = LigandSpec(
        compound_id="Mol_14056",
        pdb_id="2V5Z",
        pose_pdbqt=DATA / "poses" / "Mol_14056_2V5Z_out.pdbqt",
        net_charge=0,
    )
    ligand_mol2, ligand_frcmod = parameterize_ligand(ligand, work_dir, force=force)
    receptor_pdb = prepare_maob_receptor(work_dir)

    estimate_prefix = work_dir / "maob_mol14056_neutral_estimate"
    estimate_tleap = work_dir / "tleap_maob_estimate.in"
    write_maob_tleap(estimate_tleap, receptor_pdb, ligand_mol2, ligand_frcmod, cyf_mol2, cyf_frcmod, estimate_prefix, 0)
    run([require_tool("tleap"), "-f", str(estimate_tleap)], work_dir, "tleap_maob_estimate")
    volume_a3 = parse_box_volume(estimate_prefix.with_suffix(".inpcrd"))
    salt_pairs = salt_pairs_for_volume(volume_a3, 0.15)

    final_prefix = AMBER_OUT / "maob_mol14056"
    final_tleap = work_dir / "tleap_maob_final.in"
    write_maob_tleap(final_tleap, receptor_pdb, ligand_mol2, ligand_frcmod, cyf_mol2, cyf_frcmod, final_prefix, salt_pairs)
    run([require_tool("tleap"), "-f", str(final_tleap)], work_dir, "tleap_maob_final")

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "system": "maob_mol14056",
        "target": "MAO-B",
        "compound_id": "Mol_14056",
        "receptor": str(DATA / "receptors" / "2V5Z_chainA_fad_prepared.pdb"),
        "ligand_pose": str(ligand.pose_pdbqt),
        "ligand_net_charge": ligand.net_charge,
        "protein_force_field": "Amber ff14SB",
        "ligand_force_field": "GAFF2 with AM1-BCC charges from antechamber",
        "cofactor_parameters": {
            "residue": "CYF",
            "files": [str(cyf_mol2), str(cyf_frcmod)],
            "state": "Oxidized FAD covalently linked to Cys397 through SG-C8M thioether",
            "charge_model": "AmberTools antechamber GAFF2 AM1-BCC on capped Cys-FAD model",
            "residue_charge": -2,
        },
        "water_model": "TIP3PBOX",
        "padding_angstrom": 10.0,
        "target_salt_molar": 0.15,
        "estimated_neutral_box_volume_a3": volume_a3,
        "added_salt_pairs": salt_pairs,
        "outputs": {"prmtop": str(final_prmtop), "inpcrd": str(final_inpcrd)},
        "work_dir": str(work_dir),
    }
    (work_dir / "maob_mol14056_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return require_file(final_prmtop)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=["cyp1b1", "maob"], required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    require_tool("obabel")
    require_tool("antechamber")
    require_tool("parmchk2")
    require_tool("tleap")
    if args.target == "maob":
        prmtop = tleap_build_maob(force=args.force)
    else:
        prmtop = tleap_build_cyp1b1(force=args.force)
    print(f"Built {prmtop}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
