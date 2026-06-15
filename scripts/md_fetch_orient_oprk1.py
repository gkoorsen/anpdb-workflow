"""Fetch membrane-oriented 4DJH coordinates and orient the OPRK1 ligand pose.

The OPRK1 docking pose was generated against the original 4DJH coordinate
frame. For membrane MD, fetch the EncoMPASS/OPM-oriented 4DJH structure and
apply the same rigid-body transform to the Mol_16614 PDBQT pose.

Usage:
  python scripts/md_fetch_orient_oprk1.py
  python scripts/md_fetch_orient_oprk1.py --force
"""

from __future__ import annotations

import argparse
import json
import math
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_URL = "https://encompass.ninds.nih.gov/webfiles/database/selection/whole_structs/4djh_enc.pdb"
DEFAULT_RAW = ROOT / "data" / "md_inputs" / "receptors" / "4DJH_encompass_raw.pdb"
DEFAULT_REFERENCE = ROOT / "data" / "md_inputs" / "receptors" / "4DJH_OPRK1_clean_unoriented_reference.pdb"
FALLBACK_REFERENCE = ROOT / "output" / "docking" / "disease_link_outstanding" / "receptors" / "4DJH_OPRK1_clean.pdb"
DEFAULT_RECEPTOR_OUT = ROOT / "data" / "md_inputs" / "receptors" / "4DJH_OPRK1_opm_oriented_clean.pdb"
DEFAULT_LIGAND_IN = ROOT / "data" / "md_inputs" / "poses" / "Mol_16614_4DJH_out.pdbqt"
DEFAULT_LIGAND_OUT = ROOT / "data" / "md_inputs" / "poses" / "Mol_16614_4DJH_opm_oriented_out.pdbqt"
DEFAULT_MANIFEST = ROOT / "data" / "md_inputs" / "oprk1_orientation_manifest.json"


@dataclass(frozen=True)
class AtomRecord:
    line: str
    record: str
    serial: int
    atom_name: str
    resname: str
    chain: str
    resseq: str
    icode: str
    xyz: np.ndarray
    element: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.chain, self.resseq, self.icode, self.atom_name)


def fail(message: str) -> None:
    raise SystemExit(message)


def parse_atom_line(line: str) -> AtomRecord | None:
    if not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        serial = int(line[6:11])
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip()
        resseq = line[22:26].strip()
        icode = line[26].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
    except Exception:
        return None
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        stripped = atom_name.strip()
        if stripped and stripped[0].isdigit() and len(stripped) > 1:
            element = stripped[1]
        else:
            element = stripped[:1]
    return AtomRecord(line, line[:6].strip(), serial, atom_name, resname, chain, resseq, icode, xyz, element.upper())


def atom_field(atom_name: str, element: str) -> str:
    name = atom_name[:4]
    if len(element) == 1 and len(name) < 4:
        return f" {name:<3}"
    return f"{name:>4}"


def format_atom(record: AtomRecord, serial: int, xyz: np.ndarray) -> str:
    return (
        f"ATOM  {serial:5d} {atom_field(record.atom_name, record.element)} "
        f"{record.resname:>3} {record.chain:1}{int(record.resseq):4d}{record.icode:1}   "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          {record.element:>2}\n"
    )


def read_atoms(path: Path, chain: str | None = None, atom_name: str | None = None) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            atom = parse_atom_line(line)
            if atom is None:
                continue
            if chain is not None and atom.chain != chain:
                continue
            if atom_name is not None and atom.atom_name != atom_name:
                continue
            atoms.append(atom)
    return atoms


def download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:
        payload = response.read()
    dest.write_bytes(payload)


def kabsch_transform(reference: np.ndarray, oriented: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    if reference.shape != oriented.shape or reference.shape[0] < 10:
        fail("Need at least 10 matched atoms to compute orientation transform.")
    ref_center = reference.mean(axis=0)
    ori_center = oriented.mean(axis=0)
    ref0 = reference - ref_center
    ori0 = oriented - ori_center
    covariance = ref0.T @ ori0
    u, _s, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    translation = ori_center - ref_center @ rotation
    transformed = reference @ rotation + translation
    rmsd = math.sqrt(float(np.mean(np.sum((transformed - oriented) ** 2, axis=1))))
    return rotation, translation, rmsd


def matched_ca(reference_path: Path, oriented_path: Path, chain: str) -> tuple[np.ndarray, np.ndarray, int]:
    reference = {atom.key: atom.xyz for atom in read_atoms(reference_path, chain=chain, atom_name="CA")}
    oriented = {atom.key: atom.xyz for atom in read_atoms(oriented_path, chain=chain, atom_name="CA")}
    keys = sorted(set(reference) & set(oriented))
    if len(keys) < 10:
        fail(f"Only found {len(keys)} common CA atoms between reference and oriented 4DJH.")
    return np.array([reference[key] for key in keys]), np.array([oriented[key] for key in keys]), len(keys)


def write_clean_receptor(raw_path: Path, output_path: Path, chain: str, force: bool, source_url: str) -> int:
    if output_path.exists() and not force:
        return len(read_atoms(output_path, chain=chain))
    atoms = [atom for atom in read_atoms(raw_path, chain=chain) if atom.record == "ATOM"]
    if not atoms:
        fail(f"No ATOM records found for chain {chain} in {raw_path}.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        handle.write("REMARK  Generated by scripts/md_fetch_orient_oprk1.py\n")
        handle.write(f"REMARK  Source: {source_url}\n")
        handle.write("REMARK  Kept ATOM records for 4DJH chain A in the EncoMPASS/OPM membrane frame.\n")
        handle.write("REMARK  4DJH contains an engineered T4 lysozyme fusion segment; inspect membrane insertion.\n")
        for serial, atom in enumerate(atoms, start=1):
            handle.write(format_atom(atom, serial, atom.xyz))
        handle.write(f"TER   {len(atoms) + 1:5d}      {atoms[-1].resname:>3} {chain:1}{int(atoms[-1].resseq):4d}\n")
        handle.write("END\n")
    return len(atoms)


def transform_pdbqt(input_path: Path, output_path: Path, rotation: np.ndarray, translation: np.ndarray, force: bool) -> int:
    if output_path.exists() and not force:
        return sum(1 for line in output_path.read_text(errors="replace").splitlines() if line.startswith(("ATOM", "HETATM")))
    if not input_path.exists():
        fail(f"Missing ligand PDBQT pose: {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    transformed_atoms = 0
    out_lines: list[str] = []
    for line in input_path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
            new_xyz = xyz @ rotation + translation
            line = f"{line[:30]}{new_xyz[0]:8.3f}{new_xyz[1]:8.3f}{new_xyz[2]:8.3f}{line[54:]}"
            transformed_atoms += 1
        out_lines.append(line)
    output_path.write_text("\n".join(out_lines) + "\n")
    return transformed_atoms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--raw-output", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--receptor-output", type=Path, default=DEFAULT_RECEPTOR_OUT)
    parser.add_argument("--ligand-input", type=Path, default=DEFAULT_LIGAND_IN)
    parser.add_argument("--ligand-output", type=Path, default=DEFAULT_LIGAND_OUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    reference = args.reference
    if not reference.exists() and FALLBACK_REFERENCE.exists():
        reference = FALLBACK_REFERENCE
    if not reference.exists():
        fail(f"Missing reference 4DJH receptor: {args.reference}")

    download(args.url, args.raw_output, args.force)
    receptor_atoms = write_clean_receptor(args.raw_output, args.receptor_output, args.chain, args.force, args.url)
    reference_xyz, oriented_xyz, matched = matched_ca(reference, args.receptor_output, args.chain)
    rotation, translation, rmsd = kabsch_transform(reference_xyz, oriented_xyz)
    ligand_atoms = transform_pdbqt(args.ligand_input, args.ligand_output, rotation, translation, args.force)

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": args.url,
        "raw_output": str(args.raw_output),
        "reference": str(reference),
        "receptor_output": str(args.receptor_output),
        "ligand_input": str(args.ligand_input),
        "ligand_output": str(args.ligand_output),
        "chain": args.chain,
        "receptor_atoms": receptor_atoms,
        "transformed_ligand_atoms": ligand_atoms,
        "matched_ca_atoms": matched,
        "fit_rmsd_angstrom": rmsd,
        "rotation": rotation.tolist(),
        "translation": translation.tolist(),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

