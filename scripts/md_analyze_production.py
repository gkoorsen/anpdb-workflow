"""Analyze a completed OpenMM production MD run directory.

Outputs:
  analysis/rmsd_timeseries.csv
  analysis/rmsf_ca.csv
  analysis/contact_occupancy.csv
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
    times_ns = traj.time / 1000.0
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
    if len(lig_idx) > 0:
        lig_set = {int(idx) for idx in lig_idx}
        protein_heavy = np.array(
            [
                atom.index
                for atom in top.atoms
                if atom.residue.is_protein and atom.element is not None and atom.element.symbol != "H"
            ],
            dtype=int,
        )
        ligand_heavy = np.array(
            [
                atom.index
                for atom in top.atoms
                if atom.index in lig_set and atom.element is not None and atom.element.symbol != "H"
            ],
            dtype=int,
        )
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
    }
    (out_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
