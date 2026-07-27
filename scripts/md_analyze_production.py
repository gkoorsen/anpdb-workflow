"""Analyze a completed OpenMM production MD run directory.

Outputs:
  analysis/rmsd_timeseries.csv
  analysis/rmsf_ca.csv
  analysis/contact_occupancy.csv
  analysis/pose_retention_timeseries.csv
  analysis/ligand_geometry_timeseries.csv
  analysis/ligand_geometry_bonds.csv
  analysis/ligand_protein_hbond_occupancy.csv
  analysis/thermodynamic_timeseries.csv
  analysis/thermodynamic_qc_summary.json
  analysis/analysis_summary.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd


SOLVENT_ION_RESNAMES = {
    "HOH", "WAT", "TIP3", "SOL",
    "NA", "CL", "K", "CA", "MG", "ZN",
    "Na+", "Cl-", "K+",
}
LIPID_RESNAMES = {
    "POP", "POPC", "POPE", "POPG", "POPS", "POPA", "DPPC", "DOPC", "CHL", "CHL1",
}
COFACTOR_RESNAMES = {"HEM", "HEME", "FAD", "FMN", "NAD", "NAP"}
REPO_ROOT = Path(__file__).resolve().parents[1]
STATE_LOG_COLUMNS = {
    "Step": "step",
    "Time (ps)": "time_ps",
    "Potential Energy (kJ/mole)": "potential_energy_kJ_mol",
    "Kinetic Energy (kJ/mole)": "kinetic_energy_kJ_mol",
    "Total Energy (kJ/mole)": "total_energy_kJ_mol",
    "Temperature (K)": "temperature_K",
    "Box Volume (nm^3)": "volume_nm3",
    "Density (g/mL)": "density_g_mL",
    "Speed (ns/day)": "speed_ns_day",
}


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
        indices = np.array(
            [atom.index for atom in topology.atoms if atom.residue.name == ligand_resname],
            dtype=int,
        )
        if len(indices) > 0:
            return indices
        detected = sorted(
            {
                atom.residue.name
                for atom in topology.atoms
                if not atom.residue.is_protein
                and atom.residue.name not in SOLVENT_ION_RESNAMES | LIPID_RESNAMES | COFACTOR_RESNAMES
            }
        )
        raise SystemExit(
            f"No ligand atoms found with residue name {ligand_resname}. "
            f"Detected non-environment residue names: {detected}"
        )
    for preferred_name in ("LIG", "UNL", "UNK"):
        preferred = np.array(
            [
                atom.index
                for atom in topology.atoms
                if atom.residue.name == preferred_name
            ],
            dtype=int,
        )
        if len(preferred) > 0:
            return preferred
    excluded = SOLVENT_ION_RESNAMES | LIPID_RESNAMES | COFACTOR_RESNAMES
    indices = np.array(
        [
            atom.index
            for atom in topology.atoms
            if not atom.residue.is_protein and atom.residue.name not in excluded
        ],
        dtype=int,
    )
    if len(indices) == 0:
        raise SystemExit("No ligand atoms could be detected automatically.")
    residue_names = {
        topology.atom(int(index)).residue.name
        for index in indices
    }
    if len(residue_names) != 1:
        raise SystemExit(
            "Automatic ligand detection found multiple residue names "
            f"{sorted(residue_names)}; pass --ligand-resname explicitly."
        )
    return indices


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


def finite_or_none(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parameterized_bond_pairs(
    run_dir: Path,
    topology,
) -> tuple[np.ndarray, str]:
    system_xml = run_dir / "system.xml"
    if system_xml.exists():
        from openmm import HarmonicBondForce, XmlSerializer

        system = XmlSerializer.deserialize(system_xml.read_text())
        if system.getNumParticles() != topology.n_atoms:
            raise SystemExit(
                "system.xml particle count does not match the trajectory topology: "
                f"{system.getNumParticles()} != {topology.n_atoms}"
            )
        pairs: set[tuple[int, int]] = set()
        for force in system.getForces():
            if not isinstance(force, HarmonicBondForce):
                continue
            for index in range(force.getNumBonds()):
                atom1, atom2, _length, _k = force.getBondParameters(index)
                pairs.add(tuple(sorted((int(atom1), int(atom2)))))
        for index in range(system.getNumConstraints()):
            atom1, atom2, _distance = system.getConstraintParameters(index)
            pairs.add(tuple(sorted((int(atom1), int(atom2)))))
        if not pairs:
            raise SystemExit(f"No parameterized bonds found in {system_xml}")
        return np.array(sorted(pairs), dtype=int), "system.xml"

    config_path = run_dir / "config.toml"
    if config_path.exists():
        config = tomllib.loads(config_path.read_text())
        prmtop_value = config.get("input", {}).get("prmtop")
        if prmtop_value:
            from openmm.app import AmberPrmtopFile

            prmtop_path = Path(str(prmtop_value))
            if not prmtop_path.is_absolute():
                prmtop_path = REPO_ROOT / prmtop_path
            if not prmtop_path.exists():
                raise SystemExit(
                    f"Amber topology configured for {run_dir} is missing: {prmtop_path}"
                )
            amber_topology = AmberPrmtopFile(str(prmtop_path)).topology
            amber_atoms = list(amber_topology.atoms())
            if len(amber_atoms) != topology.n_atoms:
                raise SystemExit(
                    "Amber topology atom count does not match the trajectory topology: "
                    f"{len(amber_atoms)} != {topology.n_atoms}"
                )
            pairs = {
                tuple(sorted((atom1.index, atom2.index)))
                for atom1, atom2 in amber_topology.bonds()
            }
            if not pairs:
                raise SystemExit(f"No bonds found in Amber topology {prmtop_path}")
            return np.array(sorted(pairs), dtype=int), f"Amber prmtop: {prmtop_path}"

    raise SystemExit(
        f"No parameterized bond source found for {run_dir}. Expected system.xml "
        "or an Amber prmtop declared in config.toml."
    )


def topology_with_parameterized_bonds(topology, bond_pairs: np.ndarray):
    import mdtraj as md

    rebuilt = md.Topology()
    atom_map = {}
    for chain in topology.chains:
        new_chain = rebuilt.add_chain(chain.chain_id)
        for residue in chain.residues:
            new_residue = rebuilt.add_residue(
                residue.name,
                new_chain,
                resSeq=residue.resSeq,
                segment_id=residue.segment_id,
            )
            for atom in residue.atoms:
                new_atom = rebuilt.add_atom(
                    atom.name,
                    atom.element,
                    new_residue,
                    serial=atom.serial,
                )
                atom_map[atom.index] = new_atom
    if len(atom_map) != topology.n_atoms:
        raise SystemExit("Failed to preserve atom ordering while rebuilding topology.")
    for atom1, atom2 in bond_pairs:
        rebuilt.add_bond(atom_map[int(atom1)], atom_map[int(atom2)])
    return rebuilt


def alignment_core_indices(
    traj,
    requested_mode: str,
    custom_selection: str | None,
    membrane_core_half_thickness_nm: float,
    out_dir: Path,
    excluded_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    top = traj.topology
    excluded = {
        int(index)
        for index in (
            excluded_indices
            if excluded_indices is not None
            else np.array([], dtype=int)
        )
    }
    lipid_heavy = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.residue.name in LIPID_RESNAMES
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    membrane_center_z_nm = None
    if len(lipid_heavy) > 0:
        membrane_center_z_nm = float(traj.xyz[0, lipid_heavy, 2].mean())
    frame0_dssp = None

    if custom_selection:
        indices = np.array(
            [
                int(index)
                for index in top.select(custom_selection)
                if int(index) not in excluded
            ],
            dtype=int,
        )
        effective_mode = "custom"
        definition = custom_selection
    else:
        effective_mode = requested_mode
        if requested_mode == "auto":
            effective_mode = (
                "membrane-core-ca"
                if membrane_center_z_nm is not None
                else "global-backbone"
            )
        if effective_mode == "membrane-core-ca":
            if membrane_center_z_nm is None:
                raise SystemExit(
                    "Membrane-core alignment requested, but no lipid atoms were found."
                )
            import mdtraj as md

            frame0_dssp = md.compute_dssp(traj[0], simplified=True)[0]
            protein_ca = np.array(
                [
                    int(index)
                    for index in top.select("protein and name CA")
                    if int(index) not in excluded
                ],
                dtype=int,
            )
            relative_z = traj.xyz[0, protein_ca, 2] - membrane_center_z_nm
            indices = protein_ca[
                (np.abs(relative_z) <= membrane_core_half_thickness_nm)
                & np.array(
                    [
                        frame0_dssp[
                            top.atom(int(index)).residue.index
                        ] == "H"
                        for index in protein_ca
                    ],
                    dtype=bool,
                )
            ]
            definition = (
                "frame-0 alpha-helical protein C-alpha atoms with z coordinate "
                f"within {membrane_core_half_thickness_nm:.2f} nm of the "
                "lipid-center z"
            )
        elif effective_mode == "global-backbone":
            indices = np.array(
                [
                    int(index)
                    for index in top.select("protein and backbone")
                    if int(index) not in excluded
                ],
                dtype=int,
            )
            definition = "all protein backbone atoms"
        else:
            raise SystemExit(f"Unsupported alignment mode: {effective_mode}")

    minimum_atoms = 20 if effective_mode == "membrane-core-ca" else 3
    if len(indices) < minimum_atoms:
        raise SystemExit(
            f"Alignment core contains only {len(indices)} atoms; "
            f"at least {minimum_atoms} are required."
        )

    residue_labels = []
    rows = []
    for index in indices:
        atom = top.atom(int(index))
        residue = atom.residue
        label = f"chain{residue.chain.index}:{residue.name}{residue.resSeq}"
        residue_labels.append(label)
        rows.append(
            {
                "atom_index": int(index),
                "atom_name": atom.name,
                "residue": label,
                "frame0_secondary_structure": (
                    frame0_dssp[residue.index]
                    if frame0_dssp is not None
                    else None
                ),
                "frame0_z_nm": float(traj.xyz[0, int(index), 2]),
                "frame0_z_from_membrane_center_nm": (
                    float(traj.xyz[0, int(index), 2] - membrane_center_z_nm)
                    if membrane_center_z_nm is not None
                    else None
                ),
            }
        )
    pd.DataFrame(rows).to_csv(out_dir / "alignment_core_atoms.csv", index=False)
    metadata = {
        "requested_mode": requested_mode,
        "effective_mode": effective_mode,
        "definition": definition,
        "atom_count": int(len(indices)),
        "residue_count": int(len(set(residue_labels))),
        "membrane_center_z_nm": membrane_center_z_nm,
        "membrane_core_half_thickness_nm": (
            membrane_core_half_thickness_nm
            if effective_mode == "membrane-core-ca"
            else None
        ),
    }
    return np.asarray(indices, dtype=int), metadata


def protein_anchor_molecules(
    topology,
    excluded_indices: np.ndarray | None = None,
) -> list[set]:
    excluded = {
        int(index)
        for index in (
            excluded_indices
            if excluded_indices is not None
            else np.array([], dtype=int)
        )
    }
    protein_atom_indices = {
        int(index)
        for index in topology.select("protein")
        if int(index) not in excluded
    }
    anchors = [
        molecule
        for molecule in topology.find_molecules()
        if any(atom.index in protein_atom_indices for atom in molecule)
    ]
    if not anchors:
        raise SystemExit("No protein molecule found for PBC imaging.")
    return anchors


def image_trajectory(traj, excluded_anchor_indices: np.ndarray | None = None):
    if traj.unitcell_vectors is None:
        return traj[:]
    vectors = np.asarray(traj.unitcell_vectors)
    if vectors.size == 0 or not np.isfinite(vectors).all():
        return traj[:]
    return traj.image_molecules(
        anchor_molecules=protein_anchor_molecules(
            traj.topology,
            excluded_anchor_indices,
        ),
        inplace=False,
    )


def atom_com(traj, indices: np.ndarray) -> np.ndarray:
    masses = np.array(
        [
            float(traj.topology.atom(int(index)).element.mass)
            if traj.topology.atom(int(index)).element is not None
            else 0.0
            for index in indices
        ],
        dtype=float,
    )
    if not np.isfinite(masses).all() or masses.sum() <= 0:
        return traj.xyz[:, indices, :].mean(axis=1)
    return np.average(traj.xyz[:, indices, :], axis=1, weights=masses)


def coordinate_rmsd_A(traj, indices: np.ndarray, reference_frame: int = 0) -> np.ndarray:
    reference = traj.xyz[reference_frame, indices, :]
    deltas = traj.xyz[:, indices, :] - reference[None, :, :]
    return np.sqrt(np.mean(np.sum(deltas * deltas, axis=2), axis=1)) * 10.0


def internal_rmsd_A(traj, indices: np.ndarray) -> np.ndarray:
    import mdtraj as md

    ligand = traj.atom_slice(indices)
    return md.rmsd(ligand, ligand, frame=0) * 10.0


def analyze_ligand_geometry(
    traj,
    ligand_heavy: np.ndarray,
    parameterized_bonds: np.ndarray,
    bond_source: str,
    times_ns: np.ndarray,
    out_dir: Path,
    lower_distance_A: float = 0.8,
    upper_distance_A: float = 2.2,
) -> dict[str, object]:
    import mdtraj as md

    ligand_set = set(int(index) for index in ligand_heavy)
    bond_pairs = np.array(
        [
            [int(atom1), int(atom2)]
            for atom1, atom2 in parameterized_bonds
            if int(atom1) in ligand_set and int(atom2) in ligand_set
        ],
        dtype=int,
    )
    pdb_pairs = {
        tuple(sorted((bond.atom1.index, bond.atom2.index)))
        for bond in traj.topology.bonds
        if bond.atom1.index in ligand_set and bond.atom2.index in ligand_set
    }
    parameterized_pairs = {
        tuple(sorted((int(atom1), int(atom2))))
        for atom1, atom2 in bond_pairs
    }
    if bond_pairs.size == 0:
        return {
            "bond_source": bond_source,
            "heavy_atom_bonds": 0,
            "frames_with_implausible_heavy_atom_bond": None,
            "pdb_bonds_missing_from_parameterized_graph": len(
                pdb_pairs - parameterized_pairs
            ),
            "parameterized_bonds_missing_from_pdb_graph": len(
                parameterized_pairs - pdb_pairs
            ),
            "lower_distance_A": lower_distance_A,
            "upper_distance_A": upper_distance_A,
        }
    bond_pairs = bond_pairs.reshape((-1, 2))
    distances_A = md.compute_distances(
        traj,
        bond_pairs,
        periodic=False,
    ) * 10.0
    reference = distances_A[0]
    implausible = (distances_A < lower_distance_A) | (distances_A > upper_distance_A)
    frame_table = pd.DataFrame(
        {
            "time_ns": times_ns,
            "minimum_heavy_atom_bond_A": distances_A.min(axis=1),
            "maximum_heavy_atom_bond_A": distances_A.max(axis=1),
            "maximum_absolute_change_from_frame0_A": np.abs(
                distances_A - reference[None, :]
            ).max(axis=1),
            "implausible_heavy_atom_bond_count": implausible.sum(axis=1),
        }
    )
    frame_table.to_csv(out_dir / "ligand_geometry_timeseries.csv", index=False)

    bond_rows = []
    for column, (atom1_index, atom2_index) in enumerate(bond_pairs):
        bond_rows.append(
            {
                "atom1": str(traj.topology.atom(int(atom1_index))),
                "atom2": str(traj.topology.atom(int(atom2_index))),
                "frame0_distance_A": float(reference[column]),
                "mean_distance_A": float(distances_A[:, column].mean()),
                "minimum_distance_A": float(distances_A[:, column].min()),
                "maximum_distance_A": float(distances_A[:, column].max()),
                "implausible_frame_fraction": float(implausible[:, column].mean()),
            }
        )
    pd.DataFrame(bond_rows).to_csv(
        out_dir / "ligand_geometry_bonds.csv",
        index=False,
    )
    return {
        "bond_source": bond_source,
        "heavy_atom_bonds": int(len(bond_pairs)),
        "pdb_heavy_atom_bonds": int(len(pdb_pairs)),
        "pdb_bonds_missing_from_parameterized_graph": int(
            len(pdb_pairs - parameterized_pairs)
        ),
        "parameterized_bonds_missing_from_pdb_graph": int(
            len(parameterized_pairs - pdb_pairs)
        ),
        "frames_with_implausible_heavy_atom_bond": int(
            np.any(implausible, axis=1).sum()
        ),
        "fraction_frames_with_implausible_heavy_atom_bond": float(
            np.any(implausible, axis=1).mean()
        ),
        "minimum_observed_heavy_atom_bond_A": float(distances_A.min()),
        "maximum_observed_heavy_atom_bond_A": float(distances_A.max()),
        "maximum_absolute_change_from_frame0_A": float(
            np.abs(distances_A - reference[None, :]).max()
        ),
        "lower_distance_A": lower_distance_A,
        "upper_distance_A": upper_distance_A,
    }


def read_state_log(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "production.log"
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    frame = pd.read_csv(path)
    frame.columns = [str(column).lstrip("#").strip().strip('"') for column in frame.columns]
    available = [column for column in STATE_LOG_COLUMNS if column in frame.columns]
    out = frame[available].rename(columns=STATE_LOG_COLUMNS).copy()
    for column in out.columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    if "time_ps" in out.columns and not out["time_ps"].dropna().empty:
        out.insert(2 if "step" in out.columns else 1, "elapsed_time_ns", (out["time_ps"] - out["time_ps"].iloc[0]) / 1000.0)
    return out


def configured_barostat(run_dir: Path) -> dict[str, object]:
    path = run_dir / "config.toml"
    if not path.exists():
        return {
            "configured": None,
            "target_pressure_atm": None,
            "interval_steps": None,
        }
    config = tomllib.loads(path.read_text())
    simulation = config.get("simulation", {})
    return {
        "configured": bool(simulation.get("add_barostat", True)),
        "target_pressure_atm": float(simulation.get("pressure_atm", 1.0)),
        "interval_steps": int(simulation.get("barostat_interval_steps", 25)),
    }


def summarize_state_log(run_dir: Path, out_dir: Path) -> dict[str, object]:
    frame = read_state_log(run_dir)
    if frame.empty:
        summary = {
            "records": 0,
            "available_fields": [],
            "instantaneous_pressure_timeseries_available": False,
            "barostat": configured_barostat(run_dir),
            "metrics": {},
        }
    else:
        frame.to_csv(out_dir / "thermodynamic_timeseries.csv", index=False)
        metric_columns = [
            column
            for column in frame.columns
            if column not in {"step", "time_ps", "elapsed_time_ns", "speed_ns_day"}
        ]
        metrics: dict[str, dict[str, float | None]] = {}
        elapsed = frame.get("elapsed_time_ns")
        for column in metric_columns:
            values = frame[column].dropna()
            if values.empty:
                continue
            slope = None
            if elapsed is not None:
                paired = pd.DataFrame({"time": elapsed, "value": frame[column]}).dropna()
                if len(paired) > 1 and paired["time"].nunique() > 1:
                    slope = float(np.polyfit(paired["time"], paired["value"], 1)[0])
            metrics[column] = {
                "mean": float(values.mean()),
                "sd": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
                "min": float(values.min()),
                "max": float(values.max()),
                "linear_slope_per_ns": slope,
            }
        summary = {
            "records": int(len(frame)),
            "available_fields": list(frame.columns),
            "instantaneous_pressure_timeseries_available": False,
            "barostat": configured_barostat(run_dir),
            "metrics": metrics,
        }
    (out_dir / "thermodynamic_qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def ligand_pose_retention(traj, ligand_heavy: np.ndarray, protein_heavy: np.ndarray, lipid_indices: np.ndarray) -> pd.DataFrame:
    import mdtraj as md

    ligand_com = atom_com(traj, ligand_heavy)
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
        pocket_com = atom_com(traj, pocket_heavy)
        ligand_pocket_com_distance_A = np.linalg.norm(ligand_com - pocket_com, axis=1) * 10.0
    else:
        ligand_pocket_com_distance_A = np.full(traj.n_frames, np.nan)

    contact_data: dict[str, np.ndarray] = {}
    initial_contact_residues: set[str] = set()
    for cutoff_nm, label in [(0.4, "4A"), (0.6, "6A"), (0.8, "8A")]:
        neighbors_by_frame = md.compute_neighbors(
            traj,
            cutoff_nm,
            ligand_heavy,
            haystack_indices=protein_heavy,
            periodic=False,
        )
        residue_sets = [
            {str(traj.topology.atom(int(index)).residue) for index in frame}
            for frame in neighbors_by_frame
        ]
        contact_data[f"contact_atoms_{label}"] = np.array([len(frame) for frame in neighbors_by_frame], dtype=int)
        contact_data[f"contact_residues_{label}"] = np.array(
            [len(residues) for residues in residue_sets],
            dtype=int,
        )
        if label == "4A":
            initial_contact_residues = residue_sets[0]
            denominator = len(initial_contact_residues)
            contact_data["initial_contact_residue_fraction_4A"] = np.array(
                [
                    len(initial_contact_residues & residues) / denominator
                    if denominator
                    else np.nan
                    for residues in residue_sets
                ],
                dtype=float,
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


def ligand_protein_hbonds(
    traj,
    ligand_indices_array: np.ndarray,
    donor_acceptor_cutoff_nm: float = 0.35,
    angle_cutoff_degrees: float = 120.0,
) -> pd.DataFrame:
    import mdtraj as md

    columns = [
        "donor",
        "hydrogen",
        "acceptor",
        "ligand_role",
        "occupancy",
        "mean_donor_acceptor_distance_A_when_present",
        "donor_acceptor_cutoff_A",
        "angle_cutoff_degrees",
    ]
    ligand_set = set(int(index) for index in ligand_indices_array)
    rows = []
    try:
        hbonds = md.baker_hubbard(traj, freq=0.0, periodic=False)
    except Exception:
        hbonds = np.empty((0, 3), dtype=int)
    for donor, hydrogen, acceptor in hbonds:
        donor_atom = traj.topology.atom(int(donor))
        hydrogen_atom = traj.topology.atom(int(hydrogen))
        acceptor_atom = traj.topology.atom(int(acceptor))
        donor_lig = int(donor) in ligand_set
        acceptor_lig = int(acceptor) in ligand_set
        if donor_lig == acceptor_lig:
            continue
        distances_nm = md.compute_distances(
            traj,
            np.array([[int(donor), int(acceptor)]], dtype=int),
            periodic=False,
        )[:, 0]
        angles_rad = md.compute_angles(
            traj,
            np.array([[int(donor), int(hydrogen), int(acceptor)]], dtype=int),
            periodic=False,
        )[:, 0]
        present = (distances_nm <= donor_acceptor_cutoff_nm) & (
            angles_rad >= math.radians(angle_cutoff_degrees)
        )
        occupancy = float(present.mean())
        if occupancy <= 0:
            continue
        rows.append(
            {
                "donor": str(donor_atom),
                "hydrogen": str(hydrogen_atom),
                "acceptor": str(acceptor_atom),
                "ligand_role": "donor" if donor_lig else "acceptor",
                "occupancy": occupancy,
                "mean_donor_acceptor_distance_A_when_present": float(
                    distances_nm[present].mean() * 10.0
                ),
                "donor_acceptor_cutoff_A": donor_acceptor_cutoff_nm * 10.0,
                "angle_cutoff_degrees": angle_cutoff_degrees,
            }
        )
    return pd.DataFrame(rows, columns=columns).sort_values(
        "occupancy",
        ascending=False,
        ignore_index=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--topology", type=Path)
    parser.add_argument("--trajectory", type=Path)
    parser.add_argument("--ligand-resname")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=25,
        help=(
            "Number of trajectory frames held in memory at once. The default "
            "keeps full-resolution membrane analyses below typical 16 GB RAM."
        ),
    )
    parser.add_argument("--contact-cutoff-nm", type=float, default=0.4)
    parser.add_argument(
        "--alignment-mode",
        choices=("auto", "membrane-core-ca", "global-backbone"),
        default="auto",
        help=(
            "Auto uses a membrane-slab helical C-alpha core when lipids are present "
            "and all protein backbone atoms otherwise."
        ),
    )
    parser.add_argument(
        "--alignment-selection",
        help="Optional MDTraj selection overriding --alignment-mode.",
    )
    parser.add_argument(
        "--membrane-core-half-thickness-nm",
        type=float,
        default=1.5,
        help="Half-thickness of the membrane slab used for C-alpha core selection.",
    )
    parser.add_argument(
        "--write-imaged",
        action="store_true",
        help="Write the PBC-imaged, receptor-aligned trajectory under analysis/.",
    )
    args = parser.parse_args()

    import mdtraj as md

    run_dir = args.run_dir.resolve()
    topology_path = args.topology or run_dir / "equilibrated.pdb"
    trajectory_path = args.trajectory or run_dir / "production.dcd"
    out_dir = run_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be at least 1")

    # Load one frame to establish the authoritative atom ordering, parameterized
    # graph, imaging anchor, alignment core and all frame-0 reference values.
    raw_reference = md.load_frame(
        str(trajectory_path), 0, top=str(topology_path)
    )
    raw_ligand_indices = ligand_indices(
        raw_reference.topology, args.ligand_resname
    )
    parameterized_bonds, bond_source = parameterized_bond_pairs(
        run_dir, raw_reference.topology
    )
    parameterized_topology = topology_with_parameterized_bonds(
        raw_reference.topology,
        parameterized_bonds,
    )
    reference_unimaged = md.Trajectory(
        raw_reference.xyz.copy(),
        parameterized_topology,
        time=raw_reference.time.copy(),
        unitcell_lengths=(
            raw_reference.unitcell_lengths.copy()
            if raw_reference.unitcell_lengths is not None
            else None
        ),
        unitcell_angles=(
            raw_reference.unitcell_angles.copy()
            if raw_reference.unitcell_angles is not None
            else None
        ),
    )
    reference = image_trajectory(
        reference_unimaged,
        excluded_anchor_indices=raw_ligand_indices,
    )
    top = reference.topology
    lig_idx = ligand_indices(top, args.ligand_resname)
    lig_set = {int(index) for index in lig_idx}
    protein_backbone = np.array(
        [
            int(index)
            for index in top.select("protein and backbone")
            if int(index) not in lig_set
        ],
        dtype=int,
    )
    if len(protein_backbone) == 0:
        raise SystemExit("No protein backbone atoms found.")
    alignment_indices, alignment_metadata = alignment_core_indices(
        reference,
        args.alignment_mode,
        args.alignment_selection,
        args.membrane_core_half_thickness_nm,
        out_dir,
        excluded_indices=lig_idx,
    )

    ligand_heavy = heavy_indices(top, lig_idx)
    if len(ligand_heavy) == 0:
        raise SystemExit("The selected ligand contains no heavy atoms.")
    ca = top.select("protein and name CA")
    protein_all = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.index not in lig_set and atom.residue.is_protein
        ],
        dtype=int,
    )
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
    lipid_indices = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.residue.name in LIPID_RESNAMES
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    ligand_bond_pairs = np.array(
        [
            [int(atom1), int(atom2)]
            for atom1, atom2 in parameterized_bonds
            if int(atom1) in set(ligand_heavy)
            and int(atom2) in set(ligand_heavy)
        ],
        dtype=int,
    ).reshape((-1, 2))

    reference.superpose(reference, frame=0, atom_indices=alignment_indices)
    ref_backbone = reference.xyz[0, protein_backbone]
    ref_core = reference.xyz[0, alignment_indices]
    ref_ligand = reference.xyz[0, ligand_heavy]
    ref_ligand_com = atom_com(reference, ligand_heavy)[0]
    initial_neighbors = md.compute_neighbors(
        reference, 0.4, ligand_heavy,
        haystack_indices=protein_heavy, periodic=False,
    )[0]
    if len(initial_neighbors) == 0:
        initial_neighbors = md.compute_neighbors(
            reference, 0.8, ligand_heavy,
            haystack_indices=protein_heavy, periodic=False,
        )[0]
    pocket_heavy = np.array(sorted(set(map(int, initial_neighbors))), dtype=int)
    initial_contact_residues = {
        str(top.atom(int(index)).residue)
        for index in md.compute_neighbors(
            reference, 0.4, ligand_heavy,
            haystack_indices=protein_heavy, periodic=False,
        )[0]
    }

    with md.open(str(trajectory_path)) as handle:
        raw_frame_count = len(handle)
    selected_frame_count = (raw_frame_count + args.stride - 1) // args.stride
    log_times = read_state_log(run_dir)
    if "time_ps" in log_times.columns and len(log_times) >= raw_frame_count:
        full_times_ns = (
            log_times["time_ps"].to_numpy(dtype=float)[:: args.stride][
                :selected_frame_count
            ]
            - float(log_times["time_ps"].iloc[0])
        ) / 1000.0
    else:
        full_times_ns = np.arange(selected_frame_count, dtype=float) * 0.05 * args.stride

    rmsd_parts: list[pd.DataFrame] = []
    pose_parts: list[pd.DataFrame] = []
    geometry_parts: list[np.ndarray] = []
    contact_counts: dict[str, int] = {}
    block_labels = ("0-20", "20-40", "40-60", "60-80", "80-100")
    contact_block_counts: dict[str, dict[str, int]] = {
        label: {} for label in block_labels
    }
    ca_sum = np.zeros((len(ca), 3), dtype=np.float64)
    ca_sum_sq = np.zeros((len(ca), 3), dtype=np.float64)
    hb_accumulator: dict[tuple[str, str, str, str], list[float]] = {}
    compact_indices = np.concatenate([protein_all, lig_idx])
    compact_ligand = np.arange(len(protein_all), len(compact_indices), dtype=int)
    frame_offset = 0
    imaged_writer = None
    if args.write_imaged:
        from mdtraj.formats import DCDTrajectoryFile

        suffix = "" if args.stride == 1 else f"_stride{args.stride}"
        imaged_writer = DCDTrajectoryFile(
            str(out_dir / f"production_pbc_imaged_aligned{suffix}.dcd"), "w"
        )

    for raw_chunk in md.iterload(
        str(trajectory_path),
        top=parameterized_topology,
        chunk=args.chunk_size,
        stride=args.stride,
    ):
        chunk = image_trajectory(
            raw_chunk,
            excluded_anchor_indices=lig_idx,
        )
        chunk.superpose(reference, frame=0, atom_indices=alignment_indices)
        n_chunk = chunk.n_frames
        chunk_times = full_times_ns[frame_offset : frame_offset + n_chunk]
        if imaged_writer is not None:
            imaged_writer.write(
                chunk.xyz,
                cell_lengths=chunk.unitcell_lengths,
                cell_angles=chunk.unitcell_angles,
            )

        backbone_rmsd_chunk = np.sqrt(
            np.mean(
                np.sum(
                    (chunk.xyz[:, protein_backbone] - ref_backbone[None]) ** 2,
                    axis=2,
                ),
                axis=1,
            )
        ) * 10.0
        core_rmsd_chunk = np.sqrt(
            np.mean(
                np.sum(
                    (chunk.xyz[:, alignment_indices] - ref_core[None]) ** 2,
                    axis=2,
                ),
                axis=1,
            )
        ) * 10.0
        pose_rmsd_chunk = np.sqrt(
            np.mean(
                np.sum(
                    (chunk.xyz[:, ligand_heavy] - ref_ligand[None]) ** 2,
                    axis=2,
                ),
                axis=1,
            )
        ) * 10.0
        ligand_chunk = chunk.atom_slice(ligand_heavy)
        ligand_reference = reference.atom_slice(ligand_heavy)
        internal_chunk = md.rmsd(
            ligand_chunk, ligand_reference, frame=0
        ) * 10.0
        rmsd_parts.append(
            pd.DataFrame(
                {
                    "time_ns": chunk_times,
                    "protein_core_rmsd_A": core_rmsd_chunk,
                    "backbone_rmsd_A": backbone_rmsd_chunk,
                    "ligand_pose_rmsd_A": pose_rmsd_chunk,
                    "ligand_internal_rmsd_A": internal_chunk,
                }
            )
        )

        ca_xyz = chunk.xyz[:, ca].astype(np.float64)
        ca_sum += ca_xyz.sum(axis=0)
        ca_sum_sq += (ca_xyz * ca_xyz).sum(axis=0)
        if len(ligand_bond_pairs):
            geometry_parts.append(
                md.compute_distances(
                    chunk, ligand_bond_pairs, periodic=False
                ) * 10.0
            )

        ligand_com = atom_com(chunk, ligand_heavy)
        pocket_com = (
            atom_com(chunk, pocket_heavy)
            if len(pocket_heavy)
            else np.full((n_chunk, 3), np.nan)
        )
        pose_values: dict[str, object] = {
            "time_ns": chunk_times,
            "ligand_com_displacement_A": np.linalg.norm(
                ligand_com - ref_ligand_com[None], axis=1
            ) * 10.0,
            "ligand_pocket_com_distance_A": np.linalg.norm(
                ligand_com - pocket_com, axis=1
            ) * 10.0,
        }
        for cutoff_nm, label in ((0.4, "4A"), (0.6, "6A"), (0.8, "8A")):
            neighbors = md.compute_neighbors(
                chunk, cutoff_nm, ligand_heavy,
                haystack_indices=protein_heavy, periodic=False,
            )
            residue_sets = [
                {str(top.atom(int(index)).residue) for index in frame}
                for frame in neighbors
            ]
            pose_values[f"contact_atoms_{label}"] = [len(frame) for frame in neighbors]
            pose_values[f"contact_residues_{label}"] = [
                len(residues) for residues in residue_sets
            ]
            if label == "4A":
                denominator = len(initial_contact_residues)
                pose_values["initial_contact_residue_fraction_4A"] = [
                    (
                        len(initial_contact_residues & residues) / denominator
                        if denominator
                        else np.nan
                    )
                    for residues in residue_sets
                ]
                for residues in residue_sets:
                    for residue in residues:
                        contact_counts[residue] = contact_counts.get(residue, 0) + 1
                for local_index, residues in enumerate(residue_sets):
                    block_index = min(
                        (frame_offset + local_index) * 5
                        // selected_frame_count,
                        4,
                    )
                    block_counts = contact_block_counts[
                        block_labels[block_index]
                    ]
                    for residue in residues:
                        block_counts[residue] = block_counts.get(residue, 0) + 1
        near = md.compute_neighbors(
            chunk, 1.2, ligand_heavy,
            haystack_indices=protein_heavy, periodic=False,
        )
        minimum_distances = []
        for local_frame, neighbors in enumerate(near):
            if len(neighbors) == 0:
                minimum_distances.append(np.nan)
                continue
            delta = (
                chunk.xyz[local_frame, ligand_heavy, None, :]
                - chunk.xyz[local_frame, np.asarray(neighbors, dtype=int), :][None]
            )
            minimum_distances.append(
                float(np.sqrt(np.sum(delta * delta, axis=2)).min() * 10.0)
            )
        pose_values["protein_ligand_min_heavy_distance_A"] = minimum_distances
        if len(lipid_indices):
            pose_values["ligand_z_from_membrane_center_A"] = (
                ligand_com[:, 2] - chunk.xyz[:, lipid_indices, 2].mean(axis=1)
            ) * 10.0
        else:
            pose_values["ligand_z_from_membrane_center_A"] = np.nan

        compact = chunk.atom_slice(compact_indices)
        sasa_atom = md.shrake_rupley(compact, mode="atom")
        pose_values["ligand_sasa_nm2"] = sasa_atom[:, len(protein_all) :].sum(axis=1)
        pose_parts.append(pd.DataFrame(pose_values))

        chunk_hbonds = ligand_protein_hbonds(compact, compact_ligand)
        for row in chunk_hbonds.itertuples(index=False):
            key = (row.donor, row.hydrogen, row.acceptor, row.ligand_role)
            present = float(row.occupancy) * n_chunk
            accumulator = hb_accumulator.setdefault(key, [0.0, 0.0])
            accumulator[0] += present
            accumulator[1] += (
                present * float(row.mean_donor_acceptor_distance_A_when_present)
            )
        frame_offset += n_chunk

    if imaged_writer is not None:
        imaged_writer.close()
    if frame_offset != selected_frame_count:
        raise SystemExit(
            f"Streamed {frame_offset} frames, expected {selected_frame_count}"
        )

    rmsd_df = pd.concat(rmsd_parts, ignore_index=True)
    pose_df = pd.concat(pose_parts, ignore_index=True)
    rmsd_df.to_csv(out_dir / "rmsd_timeseries.csv", index=False)
    pose_df.to_csv(out_dir / "pose_retention_timeseries.csv", index=False)
    times_ns = rmsd_df["time_ns"].to_numpy()
    backbone_rmsd = rmsd_df["backbone_rmsd_A"].to_numpy()
    protein_core_rmsd = rmsd_df["protein_core_rmsd_A"].to_numpy()
    ligand_pose_rmsd = rmsd_df["ligand_pose_rmsd_A"].to_numpy()
    ligand_internal_rmsd = rmsd_df["ligand_internal_rmsd_A"].to_numpy()

    ca_mean = ca_sum / frame_offset
    ca_variance = np.maximum(ca_sum_sq / frame_offset - ca_mean * ca_mean, 0.0)
    rmsf = np.sqrt(ca_variance.sum(axis=1)) * 10.0
    pd.DataFrame(
        {
            "atom_index": ca,
            "residue": [str(top.atom(int(idx)).residue) for idx in ca],
            "rmsf_A": rmsf,
        }
    ).to_csv(out_dir / "rmsf_ca.csv", index=False)

    contact_df = pd.DataFrame(
        [
            {"residue": residue, "occupancy": count / frame_offset}
            for residue, count in sorted(
                contact_counts.items(), key=lambda item: (-item[1], item[0])
            )
        ]
    )
    contact_df.to_csv(out_dir / "contact_occupancy.csv", index=False)
    frame_block_index = np.minimum(
        np.arange(frame_offset, dtype=int) * 5 // frame_offset,
        4,
    )
    block_frame_counts = {
        label: int((frame_block_index == block_index).sum())
        for block_index, label in enumerate(block_labels)
    }
    contact_block_rows = []
    all_contact_residues = sorted(
        {
            residue
            for counts in contact_block_counts.values()
            for residue in counts
        }
    )
    for label in block_labels:
        denominator = block_frame_counts[label]
        for residue in all_contact_residues:
            contact_block_rows.append(
                {
                    "block_ns": label,
                    "residue": residue,
                    "occupancy": (
                        contact_block_counts[label].get(residue, 0) / denominator
                        if denominator
                        else np.nan
                    ),
                    "frames_in_block": denominator,
                }
            )
    pd.DataFrame(contact_block_rows).to_csv(
        out_dir / "contact_occupancy_20ns_blocks.csv", index=False
    )

    combined_metrics = rmsd_df.merge(pose_df, on="time_ns", how="inner")
    block_metric_columns = [
        "protein_core_rmsd_A",
        "backbone_rmsd_A",
        "ligand_pose_rmsd_A",
        "ligand_internal_rmsd_A",
        "ligand_com_displacement_A",
        "ligand_pocket_com_distance_A",
        "protein_ligand_min_heavy_distance_A",
        "contact_atoms_4A",
        "contact_residues_4A",
        "initial_contact_residue_fraction_4A",
        "ligand_sasa_nm2",
        "ligand_z_from_membrane_center_A",
    ]
    block_rows = []
    for block_index, label in enumerate(block_labels):
        block = combined_metrics.loc[frame_block_index == block_index]
        row: dict[str, object] = {
            "block_ns": label,
            "frames": int(len(block)),
        }
        for column in block_metric_columns:
            row[f"{column}_mean"] = float(block[column].mean())
            row[f"{column}_sd"] = float(block[column].std(ddof=1))
        block_rows.append(row)
    block_summary_df = pd.DataFrame(block_rows)
    block_summary_df.to_csv(
        out_dir / "metric_summary_20ns_blocks.csv", index=False
    )
    final_20 = combined_metrics.loc[frame_block_index == 4]
    final_20_summary = {
        "window_ns": "80-100",
        "frames": int(len(final_20)),
        "metrics": {
            column: {
                "mean": float(final_20[column].mean()),
                "sd": float(final_20[column].std(ddof=1)),
            }
            for column in block_metric_columns
        },
    }
    (out_dir / "final_20ns_summary.json").write_text(
        json.dumps(final_20_summary, indent=2) + "\n"
    )

    geometry_distances = (
        np.concatenate(geometry_parts, axis=0)
        if geometry_parts
        else np.empty((frame_offset, 0))
    )
    if geometry_distances.shape[1]:
        reference_bonds = geometry_distances[0]
        implausible = (geometry_distances < 0.8) | (geometry_distances > 2.2)
        pd.DataFrame(
            {
                "time_ns": times_ns,
                "minimum_heavy_atom_bond_A": geometry_distances.min(axis=1),
                "maximum_heavy_atom_bond_A": geometry_distances.max(axis=1),
                "maximum_absolute_change_from_frame0_A": np.abs(
                    geometry_distances - reference_bonds[None]
                ).max(axis=1),
                "implausible_heavy_atom_bond_count": implausible.sum(axis=1),
            }
        ).to_csv(out_dir / "ligand_geometry_timeseries.csv", index=False)
        pd.DataFrame(
            [
                {
                    "atom1": str(top.atom(int(atom1))),
                    "atom2": str(top.atom(int(atom2))),
                    "frame0_distance_A": float(reference_bonds[column]),
                    "mean_distance_A": float(geometry_distances[:, column].mean()),
                    "minimum_distance_A": float(geometry_distances[:, column].min()),
                    "maximum_distance_A": float(geometry_distances[:, column].max()),
                    "implausible_frame_fraction": float(
                        implausible[:, column].mean()
                    ),
                }
                for column, (atom1, atom2) in enumerate(ligand_bond_pairs)
            ]
        ).to_csv(out_dir / "ligand_geometry_bonds.csv", index=False)
        ligand_geometry_qc = {
            "bond_source": bond_source,
            "heavy_atom_bonds": int(len(ligand_bond_pairs)),
            "frames_with_implausible_heavy_atom_bond": int(
                np.any(implausible, axis=1).sum()
            ),
            "fraction_frames_with_implausible_heavy_atom_bond": float(
                np.any(implausible, axis=1).mean()
            ),
            "minimum_observed_heavy_atom_bond_A": float(
                geometry_distances.min()
            ),
            "maximum_observed_heavy_atom_bond_A": float(
                geometry_distances.max()
            ),
            "maximum_absolute_change_from_frame0_A": float(
                np.abs(geometry_distances - reference_bonds[None]).max()
            ),
            "lower_distance_A": 0.8,
            "upper_distance_A": 2.2,
        }
    else:
        ligand_geometry_qc = {
            "bond_source": bond_source,
            "heavy_atom_bonds": 0,
            "frames_with_implausible_heavy_atom_bond": None,
        }

    hbond_rows = []
    for (donor, hydrogen, acceptor, ligand_role), (
        present_count,
        distance_sum,
    ) in hb_accumulator.items():
        hbond_rows.append(
            {
                "donor": donor,
                "hydrogen": hydrogen,
                "acceptor": acceptor,
                "ligand_role": ligand_role,
                "occupancy": present_count / frame_offset,
                "mean_donor_acceptor_distance_A_when_present": (
                    distance_sum / present_count if present_count else np.nan
                ),
                "donor_acceptor_cutoff_A": 3.5,
                "angle_cutoff_degrees": 120.0,
            }
        )
    hbonds_df = pd.DataFrame(
        hbond_rows,
        columns=[
            "donor", "hydrogen", "acceptor", "ligand_role", "occupancy",
            "mean_donor_acceptor_distance_A_when_present",
            "donor_acceptor_cutoff_A", "angle_cutoff_degrees",
        ],
    ).sort_values("occupancy", ascending=False, ignore_index=True)
    hbonds_df.to_csv(
        out_dir / "ligand_protein_hbond_occupancy.csv", index=False
    )

    thermodynamic_qc = summarize_state_log(run_dir, out_dir)

    summary = {
        "run_dir": str(run_dir),
        "frames": frame_offset,
        "atoms": top.n_atoms,
        "chunk_size": args.chunk_size,
        "stride": args.stride,
        "time_ns_final": float(times_ns[-1]) if len(times_ns) else 0.0,
        "backbone_rmsd_A_final": float(backbone_rmsd[-1]),
        "backbone_rmsd_A_mean": float(backbone_rmsd.mean()),
        "backbone_rmsd_A_max": float(backbone_rmsd.max()),
        "protein_core_rmsd_A_final": float(protein_core_rmsd[-1]),
        "protein_core_rmsd_A_mean": float(protein_core_rmsd.mean()),
        "protein_core_rmsd_A_max": float(protein_core_rmsd.max()),
        "ligand_atoms": int(len(lig_idx)),
        "ligand_heavy_atoms": int(len(ligand_heavy)),
        "ligand_pose_rmsd_A_final": float(ligand_pose_rmsd[-1]),
        "ligand_pose_rmsd_A_mean": float(ligand_pose_rmsd.mean()),
        "ligand_pose_rmsd_A_max": float(ligand_pose_rmsd.max()),
        "ligand_internal_rmsd_A_final": float(ligand_internal_rmsd[-1]),
        "ligand_internal_rmsd_A_mean": float(ligand_internal_rmsd.mean()),
        "ligand_internal_rmsd_A_max": float(ligand_internal_rmsd.max()),
        "ligand_geometry_qc": ligand_geometry_qc,
        "pbc_imaging_applied": bool(raw_reference.unitcell_vectors is not None),
        "alignment_core": alignment_metadata,
        "contact_cutoff_nm": args.contact_cutoff_nm,
        "n_contact_residues": int(len(contact_df)),
        "ligand_com_displacement_A_final": finite_or_none(
            pose_df["ligand_com_displacement_A"].iloc[-1]
        ) if not pose_df.empty else None,
        "ligand_com_displacement_A_mean": finite_or_none(
            pose_df["ligand_com_displacement_A"].mean()
        ) if not pose_df.empty else None,
        "ligand_pocket_com_distance_A_mean": finite_or_none(
            pose_df["ligand_pocket_com_distance_A"].mean()
        ) if not pose_df.empty else None,
        "protein_ligand_min_heavy_distance_A_mean": finite_or_none(
            pose_df["protein_ligand_min_heavy_distance_A"].mean()
        ) if not pose_df.empty else None,
        "contact_atoms_4A_mean": finite_or_none(
            pose_df["contact_atoms_4A"].mean()
        ) if not pose_df.empty else None,
        "contact_atoms_6A_mean": finite_or_none(
            pose_df["contact_atoms_6A"].mean()
        ) if not pose_df.empty else None,
        "contact_atoms_8A_mean": finite_or_none(
            pose_df["contact_atoms_8A"].mean()
        ) if not pose_df.empty else None,
        "initial_contact_residue_fraction_4A_mean": (
            finite_or_none(pose_df["initial_contact_residue_fraction_4A"].mean())
            if not pose_df.empty
            else None
        ),
        "ligand_z_from_membrane_center_A_mean": finite_or_none(
            pose_df["ligand_z_from_membrane_center_A"].mean()
        ) if not pose_df.empty else None,
        "ligand_protein_hbonds_occupancy_ge_0p2": (
            int((hbonds_df["occupancy"] >= 0.2).sum()) if not hbonds_df.empty else 0
        ),
        "thermodynamic_qc": thermodynamic_qc,
    }
    (out_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
