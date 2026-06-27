"""Analyze a completed OpenMM production MD run directory.

Outputs:
  analysis/rmsd_timeseries.csv
  analysis/rmsf_ca.csv
  analysis/contact_occupancy.csv
  analysis/pose_retention_timeseries.csv
  analysis/ligand_protein_hbond_candidates.csv
  analysis/analysis_summary.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


SOLVENT_ION_RESNAMES = {
    "HOH", "WAT", "TIP3", "SOL",
    "NA", "CL", "K", "CA", "MG", "ZN",
    "Na+", "Cl-", "K+",
}
LIPID_RESNAMES = {
    "POPC", "POPE", "POPG", "POPS", "POPA", "DPPC", "DOPC", "CHL", "CHL1",
}
COFACTOR_RESNAMES = {"HEM", "HEME", "FAD", "FMN", "NAD", "NAP"}


def atom_indices_excluding(topology, excluded_resnames: set[str]) -> np.ndarray:
    return np.array(
        [
            atom.index
            for atom in topology.atoms
            if atom.residue.name not in excluded_resnames
        ],
        dtype=int,
    )


def ligand_indices(topology, ligand_resname: str | None) -> np.ndarray:
    if ligand_resname:
        return np.array(
            [atom.index for atom in topology.atoms if atom.residue.name == ligand_resname],
            dtype=int,
        )
    excluded = SOLVENT_ION_RESNAMES | LIPID_RESNAMES | COFACTOR_RESNAMES
    return np.array(
        [
            atom.index
            for atom in topology.atoms
            if not atom.residue.is_protein and atom.residue.name not in excluded
        ],
        dtype=int,
    )


def trajectory_times_ns(run_dir: Path, traj, stride: int) -> np.ndarray:
    state_log = run_dir / "production.log"
    if state_log.exists():
        log_df = pd.read_csv(state_log)
        log_df.columns = [str(column).lstrip("#").strip().strip('"') for column in log_df.columns]
        if "Time (ps)" in log_df.columns:
            times_ps = log_df["Time (ps)"].to_numpy(dtype=float)[::stride]
            times_ps = times_ps[: traj.n_frames]
            if len(times_ps) == traj.n_frames:
                return (times_ps - times_ps[0]) / 1000.0
        if "Time (ps)" in log_df.columns and len(log_df) == traj.n_frames:
            times_ps = log_df["Time (ps)"].to_numpy(dtype=float)
            return (times_ps - times_ps[0]) / 1000.0
    return traj.time / 1000.0


def heavy_indices(topology, indices: np.ndarray) -> np.ndarray:
    selected = set(int(index) for index in indices)
    return np.array(
        [
            atom.index
            for atom in topology.atoms
            if atom.index in selected and atom.element is not None and atom.element.symbol != "H"
        ],
        dtype=int,
    )


def atom_com(xyz: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return xyz[:, indices, :].mean(axis=1)


def ligand_pose_retention(traj, ligand_heavy: np.ndarray, protein_heavy: np.ndarray, lipid_indices: np.ndarray) -> pd.DataFrame:
    import mdtraj as md

    ligand_com = atom_com(traj.xyz, ligand_heavy)
    ligand_com_start = ligand_com[0]
    ligand_com_displacement_A = np.linalg.norm(ligand_com - ligand_com_start[None, :], axis=1) * 10.0

    initial_neighbors = md.compute_neighbors(
        traj[0],
        0.5,
        ligand_heavy,
        haystack_indices=protein_heavy,
        periodic=False,
    )[0]
    if len(initial_neighbors) == 0:
        initial_neighbors = md.compute_neighbors(
            traj[0],
            0.8,
            ligand_heavy,
            haystack_indices=protein_heavy,
            periodic=False,
        )[0]
    pocket_heavy = np.array(sorted(set(int(index) for index in initial_neighbors)), dtype=int)
    if len(pocket_heavy) > 0:
        pocket_com = atom_com(traj.xyz, pocket_heavy)
        ligand_pocket_com_distance_A = np.linalg.norm(ligand_com - pocket_com, axis=1) * 10.0
    else:
        ligand_pocket_com_distance_A = np.full(traj.n_frames, np.nan)

    contact_data: dict[str, np.ndarray] = {}
    for cutoff_nm, label in [(0.4, "4A"), (0.6, "6A"), (0.8, "8A")]:
        neighbors_by_frame = md.compute_neighbors(
            traj,
            cutoff_nm,
            ligand_heavy,
            haystack_indices=protein_heavy,
            periodic=False,
        )
        contact_data[f"contact_atoms_{label}"] = np.array([len(frame) for frame in neighbors_by_frame], dtype=int)
        contact_data[f"contact_residues_{label}"] = np.array(
            [len({str(traj.topology.atom(int(index)).residue) for index in frame}) for frame in neighbors_by_frame],
            dtype=int,
        )

    min_distance_A = np.empty(traj.n_frames, dtype=float)
    neighbors_12A = md.compute_neighbors(traj, 1.2, ligand_heavy, haystack_indices=protein_heavy, periodic=False)
    ligand_xyz = traj.xyz[:, ligand_heavy, :]
    for frame_index, neighbors in enumerate(neighbors_12A):
        if len(neighbors) == 0:
            min_distance_A[frame_index] = np.nan
            continue
        protein_xyz = traj.xyz[frame_index, np.array(list(neighbors), dtype=int), :]
        deltas = ligand_xyz[frame_index, :, None, :] - protein_xyz[None, :, :]
        min_distance_A[frame_index] = float(np.sqrt((deltas * deltas).sum(axis=2)).min() * 10.0)

    if len(lipid_indices) > 0:
        lipid_center_z_A = traj.xyz[:, lipid_indices, 2].mean(axis=1) * 10.0
        ligand_z_from_membrane_center_A = ligand_com[:, 2] * 10.0 - lipid_center_z_A
    else:
        lipid_center_z_A = np.full(traj.n_frames, np.nan)
        ligand_z_from_membrane_center_A = np.full(traj.n_frames, np.nan)

    data = {
        "ligand_com_displacement_A": ligand_com_displacement_A,
        "ligand_pocket_com_distance_A": ligand_pocket_com_distance_A,
        "protein_ligand_min_heavy_distance_A": min_distance_A,
        "ligand_z_from_membrane_center_A": ligand_z_from_membrane_center_A,
    }
    data.update(contact_data)
    return pd.DataFrame(data)


def ligand_protein_hbonds(traj, ligand_resname: str | None) -> pd.DataFrame:
    import mdtraj as md

    rows = []
    try:
        hbonds = md.baker_hubbard(traj, freq=0.1, periodic=False)
    except Exception:
        hbonds = np.empty((0, 3), dtype=int)
    for donor, _hydrogen, acceptor in hbonds:
        donor_atom = traj.topology.atom(int(donor))
        acceptor_atom = traj.topology.atom(int(acceptor))
        if ligand_resname:
            donor_lig = donor_atom.residue.name == ligand_resname
            acceptor_lig = acceptor_atom.residue.name == ligand_resname
        else:
            donor_lig = not donor_atom.residue.is_protein
            acceptor_lig = not acceptor_atom.residue.is_protein
        if donor_lig == acceptor_lig:
            continue
        rows.append(
            {
                "donor": str(donor_atom),
                "acceptor": str(acceptor_atom),
                "ligand_role": "donor" if donor_lig else "acceptor",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--topology", type=Path)
    parser.add_argument("--trajectory", type=Path)
    parser.add_argument("--ligand-resname")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--contact-cutoff-nm", type=float, default=0.4)
    args = parser.parse_args()

    import mdtraj as md

    run_dir = args.run_dir.resolve()
    topology_path = args.topology or run_dir / "system_solvated.pdb"
    trajectory_path = args.trajectory or run_dir / "production.dcd"
    out_dir = run_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    traj = md.load(str(trajectory_path), top=str(topology_path), stride=args.stride)
    top = traj.topology
    protein_backbone = top.select("protein and backbone")
    if len(protein_backbone) == 0:
        raise SystemExit("No protein backbone atoms found.")

    traj.superpose(traj, frame=0, atom_indices=protein_backbone)
    times_ns = trajectory_times_ns(run_dir, traj, args.stride)
    backbone_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=protein_backbone) * 10.0

    lig_idx = ligand_indices(top, args.ligand_resname)
    ligand_rmsd = None
    if len(lig_idx) > 0:
        ligand_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=lig_idx) * 10.0

    rmsd_df = pd.DataFrame({"time_ns": times_ns, "backbone_rmsd_A": backbone_rmsd})
    if ligand_rmsd is not None:
        rmsd_df["ligand_rmsd_A"] = ligand_rmsd
    rmsd_df.to_csv(out_dir / "rmsd_timeseries.csv", index=False)

    ca = top.select("protein and name CA")
    rmsf = md.rmsf(traj, traj, frame=0, atom_indices=ca) * 10.0
    rmsf_df = pd.DataFrame(
        {
            "atom_index": ca,
            "residue": [str(top.atom(int(idx)).residue) for idx in ca],
            "rmsf_A": rmsf,
        }
    )
    rmsf_df.to_csv(out_dir / "rmsf_ca.csv", index=False)

    contact_df = pd.DataFrame()
    pose_df = pd.DataFrame()
    hbonds_df = pd.DataFrame()
    if len(lig_idx) > 0:
        lig_set = {int(idx) for idx in lig_idx}
        protein_heavy = np.array(
            [
                atom.index
                for atom in top.atoms
                if atom.index not in lig_set
                and atom.residue.is_protein
                and atom.element is not None
                and atom.element.symbol != "H"
            ],
            dtype=int,
        )
        ligand_heavy = heavy_indices(top, lig_idx)
        counts: dict[str, int] = {}
        neighbors_by_frame = md.compute_neighbors(
            traj,
            args.contact_cutoff_nm,
            ligand_heavy,
            haystack_indices=protein_heavy,
        )
        for neighbors in neighbors_by_frame:
            residues = {str(top.atom(int(idx)).residue) for idx in neighbors}
            for residue in residues:
                counts[residue] = counts.get(residue, 0) + 1
        contact_df = pd.DataFrame(
            [
                {"residue": residue, "occupancy": count / traj.n_frames}
                for residue, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            ]
        )
        contact_df.to_csv(out_dir / "contact_occupancy.csv", index=False)

        lipid_indices = np.array(
            [
                atom.index
                for atom in top.atoms
                if atom.residue.name in LIPID_RESNAMES and atom.element is not None and atom.element.symbol != "H"
            ],
            dtype=int,
        )
        pose_df = ligand_pose_retention(traj, ligand_heavy, protein_heavy, lipid_indices)
        pose_df.insert(0, "time_ns", times_ns)
        pose_df.to_csv(out_dir / "pose_retention_timeseries.csv", index=False)

        compact_indices = np.array(sorted(set(protein_heavy) | set(ligand_heavy)), dtype=int)
        compact = traj.atom_slice(compact_indices)
        hbonds_df = ligand_protein_hbonds(compact, args.ligand_resname)
        hbonds_df.to_csv(out_dir / "ligand_protein_hbond_candidates.csv", index=False)

    summary = {
        "run_dir": str(run_dir),
        "frames": traj.n_frames,
        "atoms": traj.n_atoms,
        "stride": args.stride,
        "time_ns_final": float(times_ns[-1]) if len(times_ns) else 0.0,
        "backbone_rmsd_A_final": float(backbone_rmsd[-1]),
        "backbone_rmsd_A_mean": float(backbone_rmsd.mean()),
        "backbone_rmsd_A_max": float(backbone_rmsd.max()),
        "ligand_atoms": int(len(lig_idx)),
        "ligand_rmsd_A_final": float(ligand_rmsd[-1]) if ligand_rmsd is not None else None,
        "ligand_rmsd_A_mean": float(ligand_rmsd.mean()) if ligand_rmsd is not None else None,
        "ligand_rmsd_A_max": float(ligand_rmsd.max()) if ligand_rmsd is not None else None,
        "contact_cutoff_nm": args.contact_cutoff_nm,
        "n_contact_residues": int(len(contact_df)),
        "ligand_com_displacement_A_final": float(pose_df["ligand_com_displacement_A"].iloc[-1]) if not pose_df.empty else None,
        "ligand_com_displacement_A_mean": float(pose_df["ligand_com_displacement_A"].mean()) if not pose_df.empty else None,
        "ligand_pocket_com_distance_A_mean": float(pose_df["ligand_pocket_com_distance_A"].mean()) if not pose_df.empty else None,
        "protein_ligand_min_heavy_distance_A_mean": float(pose_df["protein_ligand_min_heavy_distance_A"].mean()) if not pose_df.empty else None,
        "contact_atoms_4A_mean": float(pose_df["contact_atoms_4A"].mean()) if not pose_df.empty else None,
        "contact_atoms_6A_mean": float(pose_df["contact_atoms_6A"].mean()) if not pose_df.empty else None,
        "contact_atoms_8A_mean": float(pose_df["contact_atoms_8A"].mean()) if not pose_df.empty else None,
        "ligand_z_from_membrane_center_A_mean": float(pose_df["ligand_z_from_membrane_center_A"].mean()) if not pose_df.empty else None,
        "ligand_protein_hbond_candidates_freq_ge_0p1": int(len(hbonds_df)),
    }
    (out_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
