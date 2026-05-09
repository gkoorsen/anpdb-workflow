"""Build curated Amber production input systems.

This builder is intentionally conservative. It currently creates the CYP1B1
heme system using the Shahrokh/Orendt/Yost/Cheatham P450 IC6 heme parameters
bundled with AmberTools. It does not fabricate a MAO-B covalent FAD patch; that
system needs an externally curated 8alpha-S-cysteinyl-FAD residue/parameter set.

Run from the repository root with the Docking environment:

  conda run -n Docking python scripts/md_build_amber_systems.py --target cyp1b1 --force
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "md_inputs"
AMBER_OUT = DATA / "amber_systems"
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


def explain_maob_blocker() -> None:
    message = {
        "system": "maob_mol14056",
        "ready": False,
        "blocker": "No curated covalent 8alpha-S-cysteinyl-FAD Amber residue/parameter set is available in this repo or the local AmberTools data tree.",
        "why_not_generated": [
            "Parameterizing FAD as a free ligand would give the wrong C8M valence/charges for the Cys397 thioether.",
            "Adding a manual SG-C8M bond to free-FAD GAFF parameters would not be publication-grade.",
        ],
        "required_next_input": "Curated Amber-compatible residue/library/frcmod for oxidized FAD covalently linked to Cys397, or a reviewed CHARMM-GUI/AmberTools output for the full MAO-B/FAD/Mol_14056 system.",
    }
    print(json.dumps(message, indent=2))
    raise SystemExit(2)


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
        explain_maob_blocker()

    prmtop = tleap_build_cyp1b1(force=args.force)
    print(f"Built {prmtop}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
