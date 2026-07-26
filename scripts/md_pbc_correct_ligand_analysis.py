"""PBC-correct ligand-centric compatibility analysis.

This script images each trajectory around the protein molecule, writes an
imaged DCD next to the run's analysis outputs, and computes ligand pose RMSD,
ligand SASA, and protein-ligand contact counts from the imaged trajectory.
Raw production trajectories are not modified.

New analyses should normally use ``md_analyze_production.py``, which now
performs the same PBC-aware preprocessing as part of the primary workflow.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

try:
    from md_analyze_production import (
        coordinate_rmsd_A,
        image_trajectory,
        trajectory_times_ns,
    )
except ModuleNotFoundError:
    from scripts.md_analyze_production import (
        coordinate_rmsd_A,
        image_trajectory,
        trajectory_times_ns,
    )


def ligand_indices(topology, ligand_resname: str) -> np.ndarray:
    indices = np.array([atom.index for atom in topology.atoms if atom.residue.name == ligand_resname], dtype=int)
    if len(indices) == 0:
        raise SystemExit(f"No ligand atoms found with residue name {ligand_resname}")
    return indices


def heavy_indices(topology, indices: np.ndarray) -> np.ndarray:
    index_set = set(int(idx) for idx in indices)
    return np.array(
        [
            atom.index
            for atom in topology.atoms
            if atom.index in index_set and atom.element is not None and atom.element.symbol != "H"
        ],
        dtype=int,
    )


def analyze_run(run_dir: Path, ligand_resname: str, stride: int, write_imaged: bool) -> None:
    top_path = run_dir / "equilibrated.pdb"
    traj_path = run_dir / "production.dcd"
    analysis_dir = run_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)

    traj = md.load(str(traj_path), top=str(top_path), stride=stride)
    imaged = image_trajectory(traj)

    if write_imaged:
        suffix = "" if stride == 1 else f"_stride{stride}"
        imaged.save_dcd(str(analysis_dir / f"production_pbc_imaged{suffix}.dcd"))

    top = imaged.topology
    backbone = top.select("protein and backbone")
    ligand = ligand_indices(top, ligand_resname)
    ligand_set = set(int(idx) for idx in ligand)
    protein = np.array([idx for idx in top.select("protein") if int(idx) not in ligand_set], dtype=int)
    protein_heavy = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.index not in ligand_set
            and atom.residue.is_protein
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    ligand_heavy = heavy_indices(top, ligand)

    imaged.superpose(imaged, frame=0, atom_indices=backbone)
    pose_rmsd_A = coordinate_rmsd_A(imaged, ligand_heavy)

    neighbors_4A = md.compute_neighbors(imaged, 0.4, ligand_heavy, haystack_indices=protein_heavy, periodic=False)
    neighbors_8A = md.compute_neighbors(imaged, 0.8, ligand_heavy, haystack_indices=protein_heavy, periodic=False)

    compact = imaged.atom_slice(np.concatenate([protein, ligand]))
    sasa_atom = md.shrake_rupley(compact, mode="atom")
    ligand_sasa_nm2 = sasa_atom[:, len(protein) :].sum(axis=1)

    times_ns = trajectory_times_ns(run_dir, imaged, stride)
    df = pd.DataFrame(
        {
            "time_ns": times_ns,
            "ligand_pose_rmsd_A": pose_rmsd_A,
            "ligand_sasa_nm2": ligand_sasa_nm2,
            "protein_ligand_contact_atoms_4A": [len(frame) for frame in neighbors_4A],
            "protein_ligand_contact_atoms_8A": [len(frame) for frame in neighbors_8A],
        }
    )
    df.to_csv(analysis_dir / f"pbc_corrected_ligand_metrics_stride{stride}.csv", index=False)

    summary = []
    final_time = float(times_ns[-1]) if len(times_ns) else 0.0
    windows = [(start, min(start + 20, final_time + 1e-9)) for start in np.arange(0, final_time, 20)]
    for start, end in windows:
        window = df[(df["time_ns"] >= start) & (df["time_ns"] < end)]
        if window.empty:
            continue
        summary.append(
            {
                "window_ns": f"{start:g}-{end:g}",
                "ligand_pose_rmsd_A_mean": window["ligand_pose_rmsd_A"].mean(),
                "ligand_pose_rmsd_A_max": window["ligand_pose_rmsd_A"].max(),
                "ligand_sasa_nm2_mean": window["ligand_sasa_nm2"].mean(),
                "contact_atoms_4A_mean": window["protein_ligand_contact_atoms_4A"].mean(),
                "contact_atoms_4A_min": window["protein_ligand_contact_atoms_4A"].min(),
                "contact_atoms_8A_mean": window["protein_ligand_contact_atoms_8A"].mean(),
                "contact_atoms_8A_min": window["protein_ligand_contact_atoms_8A"].min(),
            }
        )
    pd.DataFrame(summary).to_csv(analysis_dir / f"pbc_corrected_ligand_summary_stride{stride}.csv", index=False)
    print(f"Wrote PBC-corrected ligand metrics for {run_dir}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True, type=Path)
    parser.add_argument("--ligand-resname", default="LIG")
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--write-imaged", action="store_true")
    args = parser.parse_args()

    for run_dir in args.run_dir:
        analyze_run(run_dir, args.ligand_resname, args.stride, args.write_imaged)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
