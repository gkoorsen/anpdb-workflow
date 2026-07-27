"""Cluster late MD poses and prepare optional endpoint MM/PB(GB)SA inputs.

The structural analysis is performed directly from a completed OpenMM run
directory.  The final trajectory window is periodically imaged around the
receptor and aligned on the same membrane-core C-alpha definition used by
``md_analyze_production.py``.  Ligand heavy-atom pose RMSD, without fitting the
ligand to itself, is then used for hierarchical clustering.

For SGLT2 runs, this script can also convert the exact serialized OpenMM
parameterization to Amber topologies and run single-trajectory MM/PB(GB)SA.
The endpoint estimates omit configurational entropy and are intended for
within-target comparison between consistently parameterized ligands.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tomllib
import warnings

import numpy as np
import pandas as pd

from md_analyze_production import (
    LIPID_RESNAMES,
    alignment_core_indices,
    heavy_indices,
    image_trajectory,
    ligand_indices,
    parameterized_bond_pairs,
    topology_with_parameterized_bonds,
)


ALLOWED_OPENMM_FORCES = {
    "CMMotionRemover",
    "CustomExternalForce",
    "HarmonicAngleForce",
    "HarmonicBondForce",
    "MonteCarloAnisotropicBarostat",
    "MonteCarloBarostat",
    "MonteCarloMembraneBarostat",
    "NonbondedForce",
    "PeriodicTorsionForce",
}


def fail(message: str) -> None:
    raise SystemExit(message)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def load_config(run_dir: Path) -> dict:
    config_path = run_dir / "config.toml"
    if not config_path.exists():
        fail(f"Missing run configuration: {config_path}")
    return tomllib.loads(config_path.read_text())


def config_value(config: dict, *keys: str, default=None):
    value = config
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def require_run_inputs(
    run_dir: Path,
    need_solvated_topology: bool,
) -> tuple[Path, Path]:
    topology = run_dir / "equilibrated.pdb"
    trajectory = run_dir / "production.dcd"
    required = [
        topology,
        trajectory,
        run_dir / "production.log",
        run_dir / "config.toml",
        run_dir / "system.xml",
    ]
    if need_solvated_topology:
        required.append(run_dir / "system_solvated.pdb")
    missing = [path for path in required if not path.exists()]
    if missing:
        fail(
            f"Missing endpoint-analysis inputs for {run_dir}:\n"
            + "\n".join(f"  - {path}" for path in missing)
        )
    return topology, trajectory


def trajectory_frame_count(trajectory_path: Path) -> int:
    import mdtraj as md

    with md.open(str(trajectory_path)) as handle:
        return len(handle)


def trajectory_times_ns(
    run_dir: Path,
    frame_count: int,
    config: dict,
) -> np.ndarray:
    state_log = run_dir / "production.log"
    log = pd.read_csv(state_log)
    log.columns = [
        str(column).lstrip("#").strip().strip('"')
        for column in log.columns
    ]
    if "Time (ps)" in log.columns:
        numeric_times = pd.to_numeric(
            log["Time (ps)"],
            errors="coerce",
        ).dropna()
    else:
        numeric_times = pd.Series(dtype=float)
    if len(numeric_times) >= frame_count:
        times_ps = numeric_times.to_numpy(dtype=float)[:frame_count]
        return (times_ps - times_ps[0]) / 1000.0
    interval_ps = float(
        config_value(
            config,
            "output",
            "trajectory_interval_ps",
            default=50.0,
        )
    )
    return np.arange(frame_count, dtype=float) * interval_ps / 1000.0


def requested_ligand_indices(topology, ligand_resname: str | None) -> np.ndarray:
    if ligand_resname:
        requested = np.array(
            [
                atom.index
                for atom in topology.atoms
                if atom.residue.name == ligand_resname
            ],
            dtype=int,
        )
        if len(requested):
            return requested
        print(
            f"Warning: ligand residue {ligand_resname!r} was not found; "
            "falling back to unambiguous automatic detection.",
            file=sys.stderr,
        )
    return ligand_indices(topology, None)


def as_parameterized_reference(
    trajectory_path: Path,
    topology_path: Path,
    run_dir: Path,
    ligand_resname: str | None,
    out_dir: Path,
    membrane_core_half_thickness_nm: float,
):
    import mdtraj as md

    raw_reference = md.load_frame(
        str(trajectory_path),
        0,
        top=str(topology_path),
    )
    raw_ligand = requested_ligand_indices(
        raw_reference.topology,
        ligand_resname,
    )
    parameterized_bonds, bond_source = parameterized_bond_pairs(
        run_dir,
        raw_reference.topology,
    )
    parameterized_topology = topology_with_parameterized_bonds(
        raw_reference.topology,
        parameterized_bonds,
    )
    reference = md.Trajectory(
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
        reference,
        excluded_anchor_indices=raw_ligand,
    )
    ligand = requested_ligand_indices(reference.topology, ligand_resname)
    alignment, alignment_metadata = alignment_core_indices(
        reference,
        "auto",
        None,
        membrane_core_half_thickness_nm,
        out_dir,
        excluded_indices=ligand,
    )
    reference.superpose(reference, frame=0, atom_indices=alignment)
    ligand_heavy = heavy_indices(reference.topology, ligand)
    if not len(ligand_heavy):
        fail("The selected ligand has no heavy atoms.")

    ligand_set = set(map(int, ligand))
    protein = np.array(
        [
            atom.index
            for atom in reference.topology.atoms
            if atom.index not in ligand_set and atom.residue.is_protein
        ],
        dtype=int,
    )
    if not len(protein):
        fail("No receptor protein atoms were found.")
    complex_indices = np.concatenate([protein, ligand])
    if not np.array_equal(complex_indices, np.sort(complex_indices)):
        fail(
            "The receptor and ligand are not ordered as receptor-then-ligand "
            "in the trajectory. This topology requires manual review before "
            "MM/PB(GB)SA."
        )
    lipid_heavy = np.array(
        [
            atom.index
            for atom in reference.topology.atoms
            if atom.residue.name in LIPID_RESNAMES
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    return {
        "reference": reference,
        "topology": parameterized_topology,
        "parameterized_bonds": parameterized_bonds,
        "bond_source": bond_source,
        "ligand": ligand,
        "ligand_heavy": ligand_heavy,
        "protein": protein,
        "complex": complex_indices,
        "lipid_heavy": lipid_heavy,
        "alignment": alignment,
        "alignment_metadata": alignment_metadata,
    }


def late_frame_indices(
    times_ns: np.ndarray,
    late_window_ns: float,
    stride: int,
) -> np.ndarray:
    if late_window_ns <= 0:
        fail("--late-window-ns must be positive.")
    if stride < 1:
        fail("--cluster-stride must be at least 1.")
    start_ns = float(times_ns[-1]) - late_window_ns
    selected = np.flatnonzero(times_ns >= start_ns - 1e-9)[::stride]
    if len(selected) < 2:
        fail(
            f"Only {len(selected)} frame(s) were selected from the final "
            f"{late_window_ns:g} ns."
        )
    return selected


def compact_topology_with_source_numbering(
    run_dir: Path,
    topology,
    complex_indices: np.ndarray,
    protein_indices: np.ndarray,
    parameterized_bonds: np.ndarray,
    out_dir: Path,
):
    import mdtraj as md

    selected = set(map(int, complex_indices))
    protein_set = set(map(int, protein_indices))
    system_protein_residues = []
    for residue in topology.residues:
        if any(atom.index in protein_set for atom in residue.atoms):
            system_protein_residues.append(residue)

    source_path = run_dir / "receptor_fixed.pdb"
    source_residues = []
    if source_path.exists():
        source = md.load_frame(str(source_path), 0)
        source_residues = [
            residue
            for residue in source.topology.residues
            if residue.is_protein
        ]

    source_mapping_valid = (
        len(source_residues) == len(system_protein_residues)
        and all(
            source.name == system.name
            for source, system in zip(
                source_residues,
                system_protein_residues,
            )
        )
    )

    rebuilt = md.Topology()
    atom_map = {}
    chain_map = {}
    residue_map_rows = []
    protein_ordinal = 0
    for residue in topology.residues:
        residue_atoms = [
            atom
            for atom in residue.atoms
            if atom.index in selected
        ]
        if not residue_atoms:
            continue
        is_protein = any(atom.index in protein_set for atom in residue_atoms)
        source_residue = None
        if is_protein:
            if source_mapping_valid:
                source_residue = source_residues[protein_ordinal]
            protein_ordinal += 1
        source_chain = source_residue.chain if source_residue else residue.chain
        chain_key = (
            "protein" if is_protein else "ligand",
            source_chain.index,
        )
        if chain_key not in chain_map:
            chain_map[chain_key] = rebuilt.add_chain(source_chain.chain_id)
        new_residue = rebuilt.add_residue(
            residue.name,
            chain_map[chain_key],
            resSeq=(
                source_residue.resSeq
                if source_residue is not None
                else residue.resSeq
            ),
            segment_id=residue.segment_id,
        )
        residue_map_rows.append(
            {
                "system_residue_index": residue.index,
                "system_residue": f"{residue.name}{residue.resSeq}",
                "source_chain": (
                    source_residue.chain.chain_id
                    if source_residue is not None
                    else residue.chain.chain_id
                ),
                "source_residue_number": (
                    source_residue.resSeq
                    if source_residue is not None
                    else residue.resSeq
                ),
                "source_residue": (
                    f"{source_residue.name}{source_residue.resSeq}"
                    if source_residue is not None
                    else f"{residue.name}{residue.resSeq}"
                ),
                "mapping_source": (
                    "receptor_fixed.pdb"
                    if source_residue is not None
                    else "trajectory_topology"
                ),
            }
        )
        for atom in residue_atoms:
            atom_map[atom.index] = rebuilt.add_atom(
                atom.name,
                atom.element,
                new_residue,
                serial=atom.serial,
            )

    for atom1, atom2 in parameterized_bonds:
        atom1 = int(atom1)
        atom2 = int(atom2)
        if atom1 in atom_map and atom2 in atom_map:
            rebuilt.add_bond(atom_map[atom1], atom_map[atom2])

    if len(atom_map) != len(complex_indices):
        fail("Failed to preserve complex atom ordering in the compact topology.")
    pd.DataFrame(residue_map_rows).to_csv(
        out_dir / "residue_number_map.csv",
        index=False,
    )
    if source_residues and not source_mapping_valid:
        print(
            "Warning: receptor_fixed.pdb could not be mapped one-to-one to "
            "the system protein; medoid PDBs retain trajectory numbering.",
            file=sys.stderr,
        )
    return rebuilt, source_mapping_valid


def collect_selected_frames(
    trajectory_path: Path,
    selected_indices: np.ndarray,
    prepared: dict,
    chunk_size: int,
):
    import mdtraj as md

    if chunk_size < 1:
        fail("--chunk-size must be at least 1.")
    selected_indices = np.asarray(selected_indices, dtype=int)
    selected_set = set(map(int, selected_indices))
    frame_order = {int(frame): index for index, frame in enumerate(selected_indices)}
    ligand_xyz = [None] * len(selected_indices)
    complex_xyz = [None] * len(selected_indices)
    membrane_center_z = np.full(len(selected_indices), np.nan, dtype=float)

    offset = 0
    for raw_chunk in md.iterload(
        str(trajectory_path),
        top=prepared["topology"],
        chunk=chunk_size,
    ):
        global_indices = [
            index
            for index in range(offset, offset + raw_chunk.n_frames)
            if index in selected_set
        ]
        if global_indices:
            local_indices = np.array(global_indices, dtype=int) - offset
            chunk = raw_chunk[local_indices]
            chunk = image_trajectory(
                chunk,
                excluded_anchor_indices=prepared["ligand"],
            )
            chunk.superpose(
                prepared["reference"],
                frame=0,
                atom_indices=prepared["alignment"],
            )
            for local, global_index in enumerate(global_indices):
                destination = frame_order[global_index]
                ligand_xyz[destination] = chunk.xyz[
                    local,
                    prepared["ligand_heavy"],
                    :,
                ].copy()
                complex_xyz[destination] = chunk.xyz[
                    local,
                    prepared["complex"],
                    :,
                ].copy()
                if len(prepared["lipid_heavy"]):
                    membrane_center_z[destination] = float(
                        chunk.xyz[
                            local,
                            prepared["lipid_heavy"],
                            2,
                        ].mean()
                    )
        offset += raw_chunk.n_frames

    if any(value is None for value in ligand_xyz + complex_xyz):
        fail("Not all requested frames were recovered from the trajectory.")
    return (
        np.asarray(ligand_xyz, dtype=np.float32),
        np.asarray(complex_xyz, dtype=np.float32),
        membrane_center_z,
    )


def renumber_clusters(labels: np.ndarray) -> np.ndarray:
    counts = {
        int(label): int((labels == label).sum())
        for label in np.unique(labels)
    }
    ordered = sorted(
        counts,
        key=lambda label: (
            -counts[label],
            int(np.flatnonzero(labels == label)[0]),
        ),
    )
    mapping = {old: new for new, old in enumerate(ordered, start=1)}
    return np.array([mapping[int(label)] for label in labels], dtype=int)


def silhouette_from_distance_matrix(
    distance_matrix: np.ndarray,
    labels: np.ndarray,
) -> float | None:
    unique = np.unique(labels)
    if len(unique) < 2 or len(unique) >= len(labels):
        return None
    values = []
    for index, label in enumerate(labels):
        same = np.flatnonzero(labels == label)
        same = same[same != index]
        if not len(same):
            values.append(0.0)
            continue
        a_value = float(distance_matrix[index, same].mean())
        b_value = min(
            float(
                distance_matrix[
                    index,
                    np.flatnonzero(labels == other),
                ].mean()
            )
            for other in unique
            if other != label
        )
        denominator = max(a_value, b_value)
        values.append(
            (b_value - a_value) / denominator
            if denominator > 0
            else 0.0
        )
    return float(np.mean(values))


def cluster_ligand_poses(
    ligand_xyz_nm: np.ndarray,
    frame_indices: np.ndarray,
    times_ns: np.ndarray,
    distance_threshold_A: float,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, np.ndarray]:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist, squareform

    if distance_threshold_A <= 0:
        fail("--cluster-threshold-A must be positive.")
    flattened = ligand_xyz_nm.reshape((len(ligand_xyz_nm), -1))
    condensed_A = (
        pdist(flattened, metric="euclidean")
        / math.sqrt(ligand_xyz_nm.shape[1])
        * 10.0
    )
    linkage_matrix = linkage(condensed_A, method="average")
    labels = renumber_clusters(
        fcluster(
            linkage_matrix,
            t=distance_threshold_A,
            criterion="distance",
        )
    )
    distance_matrix = squareform(condensed_A)
    medoid_by_cluster = {}
    rows = []
    for cluster_id in sorted(np.unique(labels)):
        members = np.flatnonzero(labels == cluster_id)
        within = distance_matrix[np.ix_(members, members)]
        medoid_local = int(np.argmin(within.sum(axis=1)))
        medoid = int(members[medoid_local])
        medoid_by_cluster[int(cluster_id)] = medoid
        member_distances = distance_matrix[medoid, members]
        rows.append(
            {
                "cluster_id": int(cluster_id),
                "frames": int(len(members)),
                "fraction": float(len(members) / len(labels)),
                "medoid_selected_index": medoid,
                "medoid_frame_index": int(frame_indices[medoid]),
                "medoid_time_ns": float(times_ns[medoid]),
                "mean_distance_to_medoid_A": float(member_distances.mean()),
                "median_distance_to_medoid_A": float(
                    np.median(member_distances)
                ),
                "maximum_distance_to_medoid_A": float(
                    member_distances.max()
                ),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["frames", "cluster_id"],
        ascending=[False, True],
    )
    assignment_rows = []
    for index, cluster_id in enumerate(labels):
        medoid = medoid_by_cluster[int(cluster_id)]
        assignment_rows.append(
            {
                "selected_index": index,
                "frame_index": int(frame_indices[index]),
                "time_ns": float(times_ns[index]),
                "cluster_id": int(cluster_id),
                "distance_to_cluster_medoid_A": float(
                    distance_matrix[index, medoid]
                ),
                "is_medoid": bool(index == medoid),
            }
        )
    assignments = pd.DataFrame(assignment_rows)
    assignments.to_csv(out_dir / "cluster_assignments.csv", index=False)
    summary.to_csv(out_dir / "cluster_summary.csv", index=False)
    np.save(out_dir / "cluster_linkage.npy", linkage_matrix)

    dominant = summary.iloc[0]
    manifest = {
        "method": "average-linkage hierarchical clustering",
        "distance_metric": (
            "ligand heavy-atom pose RMSD after receptor-core alignment; "
            "the ligand was not fitted to itself"
        ),
        "distance_threshold_A": distance_threshold_A,
        "selected_frames": int(len(labels)),
        "clusters": int(len(summary)),
        "dominant_cluster_id": int(dominant["cluster_id"]),
        "dominant_cluster_fraction": float(dominant["fraction"]),
        "dominant_cluster_medoid_frame_index": int(
            dominant["medoid_frame_index"]
        ),
        "dominant_cluster_medoid_time_ns": float(
            dominant["medoid_time_ns"]
        ),
        "mean_silhouette": silhouette_from_distance_matrix(
            distance_matrix,
            labels,
        ),
    }
    write_json(out_dir / "clustering_manifest.json", manifest)
    return assignments, summary, manifest, labels


def write_medoid_structures(
    compact_xyz_nm: np.ndarray,
    compact_topology,
    summary: pd.DataFrame,
    out_dir: Path,
    maximum_structures: int,
) -> None:
    import mdtraj as md

    if maximum_structures < 1:
        fail("--max-medoid-structures must be at least 1.")
    medoid_dir = out_dir / "cluster_medoids"
    medoid_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for row in summary.head(maximum_structures).itertuples(index=False):
        index = int(row.medoid_selected_index)
        path = medoid_dir / f"cluster_{int(row.cluster_id)}_medoid.pdb"
        md.Trajectory(
            compact_xyz_nm[index : index + 1],
            compact_topology,
        ).save_pdb(str(path))
        written.append(path)
    shutil.copyfile(
        written[0],
        out_dir / "dominant_cluster_medoid.pdb",
    )


def validate_openmm_system_for_conversion(system) -> dict:
    force_rows = []
    for force in system.getForces():
        name = type(force).__name__
        if name not in ALLOWED_OPENMM_FORCES:
            fail(
                f"OpenMM-to-Amber conversion has not been validated for {name}."
            )
        row = {"force": name}
        if name == "CustomExternalForce":
            parameters = {}
            for index in range(force.getNumGlobalParameters()):
                parameter_name = force.getGlobalParameterName(index)
                parameters[parameter_name] = float(
                    force.getGlobalParameterDefaultValue(index)
                )
            row["global_parameters"] = parameters
            nonzero = {
                key: value
                for key, value in parameters.items()
                if abs(value) > 1e-12
            }
            if nonzero:
                fail(
                    "The serialized CustomExternalForce has nonzero default "
                    f"parameters and cannot be omitted safely: {nonzero}"
                )
        force_rows.append(row)
    nonbonded = [
        force
        for force in system.getForces()
        if type(force).__name__ == "NonbondedForce"
    ]
    if len(nonbonded) != 1:
        fail(
            "Expected exactly one standard OpenMM NonbondedForce for "
            f"conversion; found {len(nonbonded)}."
        )
    return {
        "forces": force_rows,
        "particles": int(system.getNumParticles()),
        "constraints": int(system.getNumConstraints()),
    }


def assign_amber_radii(structure, radii_set: str) -> None:
    from parmed.tools.actions import changeRadii

    changeRadii(structure, radii_set).execute()
    invalid = [
        atom.idx
        for atom in structure.atoms
        if atom.solvent_radius <= 0 or atom.screen <= 0
    ]
    if invalid:
        fail(
            f"Amber {radii_set} radii were not assigned to "
            f"{len(invalid)} atom(s)."
        )


def subset_structure(structure, indices: np.ndarray):
    mask = np.zeros(len(structure.atoms), dtype=bool)
    mask[np.asarray(indices, dtype=int)] = True
    subset = structure[mask]
    subset.box = None
    return subset


def prepare_amber_topologies(
    run_dir: Path,
    prepared: dict,
    out_dir: Path,
    radii_set: str,
) -> dict:
    import parmed as pmd
    from openmm import XmlSerializer
    from openmm.app import PDBFile

    system_path = run_dir / "system.xml"
    topology_path = run_dir / "system_solvated.pdb"
    system = XmlSerializer.deserialize(system_path.read_text())
    force_manifest = validate_openmm_system_for_conversion(system)
    pdb = PDBFile(str(topology_path))
    if system.getNumParticles() != pdb.topology.getNumAtoms():
        fail(
            "system.xml and system_solvated.pdb atom counts differ: "
            f"{system.getNumParticles()} != {pdb.topology.getNumAtoms()}"
        )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        full = pmd.openmm.load_topology(
            pdb.topology,
            system,
            xyz=pdb.positions,
        )
    warning_messages = [str(item.message) for item in caught]
    unsupported = [
        message
        for message in warning_messages
        if "not supported" in message.lower()
        or "unknown functional" in message.lower()
    ]
    if unsupported:
        fail(
            "ParmEd reported unsupported OpenMM terms:\n"
            + "\n".join(f"  - {message}" for message in unsupported)
        )
    if len(full.atoms) != system.getNumParticles():
        fail("ParmEd changed the atom count during OpenMM conversion.")

    complex_structure = subset_structure(full, prepared["complex"])
    receptor_structure = subset_structure(full, prepared["protein"])
    ligand_structure = subset_structure(full, prepared["ligand"])
    if len(complex_structure.atoms) != (
        len(receptor_structure.atoms) + len(ligand_structure.atoms)
    ):
        fail("Complex, receptor and ligand topology atom counts are inconsistent.")

    for structure in (
        complex_structure,
        receptor_structure,
        ligand_structure,
    ):
        assign_amber_radii(structure, radii_set)

    topology_dir = out_dir / "mmpbsa"
    topology_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "complex": topology_dir / "complex.prmtop",
        "receptor": topology_dir / "receptor.prmtop",
        "ligand": topology_dir / "ligand.prmtop",
    }
    complex_structure.save(str(paths["complex"]), overwrite=True)
    receptor_structure.save(str(paths["receptor"]), overwrite=True)
    ligand_structure.save(str(paths["ligand"]), overwrite=True)

    reloaded_counts = {}
    for name, path in paths.items():
        loaded = pmd.load_file(str(path))
        reloaded_counts[name] = len(loaded.atoms)
        expected = len(
            {
                "complex": complex_structure,
                "receptor": receptor_structure,
                "ligand": ligand_structure,
            }[name].atoms
        )
        if reloaded_counts[name] != expected:
            fail(
                f"Reloaded {name} prmtop has {reloaded_counts[name]} atoms; "
                f"expected {expected}."
            )

    charges = {
        "complex": float(sum(atom.charge for atom in complex_structure.atoms)),
        "receptor": float(sum(atom.charge for atom in receptor_structure.atoms)),
        "ligand": float(sum(atom.charge for atom in ligand_structure.atoms)),
    }
    if abs(charges["complex"] - charges["receptor"] - charges["ligand"]) > 1e-5:
        fail("Complex charge does not equal receptor plus ligand charge.")

    manifest = {
        "source_system": str(system_path),
        "source_topology": str(topology_path),
        "conversion": "parmed.openmm.load_topology",
        "radii_set": radii_set,
        "atom_counts": reloaded_counts,
        "charges_e": charges,
        "parmed_warnings": warning_messages,
        "openmm_system": force_manifest,
        "omitted_openmm_terms": [
            "Monte Carlo barostat: no potential-energy contribution",
            "CMMotionRemover: no potential-energy contribution",
            (
                "CustomExternalForce equilibration restraint: accepted only "
                "after confirming all default global parameters are zero"
            ),
        ],
    }
    write_json(topology_dir / "topology_conversion_manifest.json", manifest)
    return {
        "directory": topology_dir,
        "paths": paths,
        "manifest": manifest,
    }


def evenly_spaced_indices(length: int, count: int) -> np.ndarray:
    if count < 1:
        fail("--mmpbsa-snapshots must be at least 1.")
    count = min(length, count)
    return np.unique(
        np.rint(np.linspace(0, length - 1, count)).astype(int)
    )


def write_mmpbsa_trajectory(
    compact_xyz_nm: np.ndarray,
    compact_topology,
    membrane_center_z_nm: np.ndarray,
    selected_frame_indices: np.ndarray,
    selected_times_ns: np.ndarray,
    out_dir: Path,
    snapshot_count: int,
) -> tuple[Path, pd.DataFrame]:
    import mdtraj as md

    selected = evenly_spaced_indices(len(compact_xyz_nm), snapshot_count)
    coordinates = compact_xyz_nm[selected].copy()
    centers = membrane_center_z_nm[selected]
    if not np.isfinite(centers).all():
        fail(
            "Implicit-membrane MM/PBSA requires lipid coordinates so that "
            "the membrane center can be translated to z=0."
        )
    coordinates[:, :, 2] -= centers[:, None]
    trajectory = md.Trajectory(
        coordinates,
        compact_topology,
        time=selected_times_ns[selected] * 1000.0,
    )
    mmpbsa_dir = out_dir / "mmpbsa"
    mmpbsa_dir.mkdir(parents=True, exist_ok=True)
    trajectory_path = mmpbsa_dir / "complex_snapshots.nc"
    trajectory.save_netcdf(str(trajectory_path))
    trajectory[0].save_pdb(str(mmpbsa_dir / "complex_snapshot_001.pdb"))
    frames = pd.DataFrame(
        {
            "snapshot": np.arange(1, len(selected) + 1),
            "selected_index": selected,
            "trajectory_frame_index": selected_frame_indices[selected],
            "time_ns": selected_times_ns[selected],
            "original_membrane_center_z_nm": centers,
            "translated_membrane_center_z_nm": 0.0,
        }
    )
    frames.to_csv(mmpbsa_dir / "snapshot_frames.csv", index=False)
    return trajectory_path, frames


def write_mmpbsa_input(
    path: Path,
    frame_count: int,
    salt_molar: float,
    solute_dielectric: float,
    membrane_dielectric: float,
    membrane_thickness_A: float,
) -> None:
    if not (0 < solute_dielectric <= membrane_dielectric < 80.0):
        fail(
            "Require 0 < solute dielectric <= membrane dielectric < 80."
        )
    path.write_text(
        "\n".join(
            [
                "Single-trajectory endpoint MM/PB(GB)SA for SGLT2",
                "&general",
                f"  startframe=1, endframe={frame_count}, interval=1,",
                "  keep_files=0, use_sander=1, verbose=1,",
                "/",
                "&gb",
                f"  igb=5, saltcon={salt_molar:.4f},",
                "/",
                "&pb",
                "  radiopt=0,",
                f"  indi={solute_dielectric:.4f}, exdi=80.0,",
                f"  istrng={salt_molar:.4f},",
                "  fillratio=1.25, inp=2,",
                "  sasopt=0, solvopt=2, ipb=1, bcopt=10,",
                "  nfocus=1, linit=1000, eneopt=1,",
                "  cutfd=7.0, cutnb=99.0, maxarcdot=15000, npbverb=1,",
                "  memopt=1,",
                f"  emem={membrane_dielectric:.4f},",
                "  mctrdz=0.0,",
                f"  mthick={membrane_thickness_A:.4f}, poretype=1,",
                "/",
                "",
            ]
        )
    )


def mmpbsa_command(
    mmpbsa_dir: Path,
    mpi_ranks: int,
    executable: str | None,
) -> list[str]:
    serial = executable or shutil.which("MMPBSA.py")
    mpi = shutil.which("MMPBSA.py.MPI")
    if mpi_ranks < 1:
        fail("--mpi-ranks must be at least 1.")
    base = [
        "-O",
        "-i",
        "mmpbsa.in",
        "-o",
        "FINAL_RESULTS_MMPBSA.dat",
        "-eo",
        "FINAL_RESULTS_MMPBSA.csv",
        "-sp",
        "complex.prmtop",
        "-cp",
        "complex.prmtop",
        "-rp",
        "receptor.prmtop",
        "-lp",
        "ligand.prmtop",
        "-y",
        "complex_snapshots.nc",
    ]
    if mpi_ranks == 1:
        if not serial:
            fail(
                "MMPBSA.py was not found. Create/activate "
                "environment-md-endpoint.yml or use --mmpbsa-executable."
            )
        return [serial, *base]
    launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if not launcher or not mpi:
        fail(
            "Parallel execution requires mpirun/mpiexec and MMPBSA.py.MPI."
        )
    command = [launcher]
    if os.geteuid() == 0 and Path(launcher).name == "mpirun":
        command.append("--allow-run-as-root")
    return [*command, "-np", str(mpi_ranks), mpi, *base]


def parse_mmpbsa_results(path: Path) -> list[dict]:
    if not path.exists():
        return []
    model = None
    in_difference = False
    rows = []
    binding_pattern = re.compile(
        r"DELTA G binding\s*=\s*"
        r"([-+0-9.Ee]+)\s*\+/-\s*"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)"
    )
    total_pattern = re.compile(
        r"^\s*DELTA TOTAL\s+"
        r"([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)"
    )
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("GENERALIZED BORN"):
            model = "GB"
            in_difference = False
        elif stripped.startswith("POISSON BOLTZMANN"):
            model = "PB-implicit-membrane"
            in_difference = False
        elif stripped.startswith("Differences"):
            in_difference = True
        if not model:
            continue
        match = binding_pattern.search(line)
        if match:
            rows.append(
                {
                    "model": model,
                    "delta_endpoint_kcal_mol": float(match.group(1)),
                    "snapshot_sd_kcal_mol": float(match.group(2)),
                    "snapshot_sem_kcal_mol": float(match.group(3)),
                }
            )
            in_difference = False
            continue
        match = total_pattern.match(line) if in_difference else None
        if match:
            rows.append(
                {
                    "model": model,
                    "delta_endpoint_kcal_mol": float(match.group(1)),
                    "snapshot_sd_kcal_mol": float(match.group(2)),
                    "snapshot_sem_kcal_mol": float(match.group(3)),
                }
            )
            in_difference = False
    deduplicated = {}
    for row in rows:
        deduplicated[row["model"]] = row
    return list(deduplicated.values())


def run_mmpbsa(
    mmpbsa_dir: Path,
    mpi_ranks: int,
    executable: str | None,
) -> list[dict]:
    command = mmpbsa_command(mmpbsa_dir, mpi_ranks, executable)
    write_json(
        mmpbsa_dir / "mmpbsa_command.json",
        {"argv": command, "working_directory": str(mmpbsa_dir)},
    )
    with (mmpbsa_dir / "mmpbsa.stdout.log").open("w") as log:
        subprocess.run(
            command,
            cwd=mmpbsa_dir,
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    rows = parse_mmpbsa_results(
        mmpbsa_dir / "FINAL_RESULTS_MMPBSA.dat"
    )
    if not rows:
        fail(
            "MMPBSA.py completed but no endpoint estimate could be parsed from "
            f"{mmpbsa_dir / 'FINAL_RESULTS_MMPBSA.dat'}"
        )
    pd.DataFrame(rows).to_csv(
        mmpbsa_dir / "mmpbsa_summary.csv",
        index=False,
    )
    write_json(
        mmpbsa_dir / "mmpbsa_summary.json",
        {
            "models": rows,
            "interpretation": (
                "Single-trajectory endpoint estimates without configurational "
                "entropy; compare only consistently parameterized ligands at "
                "the same target."
            ),
        },
    )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--ligand-resname")
    parser.add_argument("--late-window-ns", type=float, default=20.0)
    parser.add_argument("--cluster-stride", type=int, default=1)
    parser.add_argument("--cluster-threshold-A", type=float, default=2.0)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--max-medoid-structures", type=int, default=5)
    parser.add_argument(
        "--skip-clustering",
        action="store_true",
        help="Prepare/run endpoint energies without rewriting cluster outputs.",
    )
    parser.add_argument(
        "--prepare-mmpbsa",
        action="store_true",
        help="Prepare Amber topologies and late-stage snapshots.",
    )
    parser.add_argument(
        "--run-mmpbsa",
        action="store_true",
        help="Prepare inputs and execute MMPBSA.py.",
    )
    parser.add_argument("--mmpbsa-snapshots", type=int, default=100)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mmpbsa-executable")
    parser.add_argument("--amber-radii", default="mbondi2")
    parser.add_argument("--salt-molar", type=float, default=0.15)
    parser.add_argument("--solute-dielectric", type=float, default=4.0)
    parser.add_argument("--membrane-dielectric", type=float, default=7.0)
    parser.add_argument("--membrane-thickness-A", type=float, default=40.0)
    parser.add_argument(
        "--membrane-core-half-thickness-nm",
        type=float,
        default=1.5,
    )
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    prepare_energy = args.prepare_mmpbsa or args.run_mmpbsa
    topology_path, trajectory_path = require_run_inputs(
        run_dir,
        need_solvated_topology=prepare_energy,
    )
    config = load_config(run_dir)
    out_dir = run_dir / "endpoint_analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    frame_count = trajectory_frame_count(trajectory_path)
    times_ns = trajectory_times_ns(run_dir, frame_count, config)
    selected_frames = late_frame_indices(
        times_ns,
        args.late_window_ns,
        args.cluster_stride,
    )
    selected_times = times_ns[selected_frames]
    prepared = as_parameterized_reference(
        trajectory_path,
        topology_path,
        run_dir,
        args.ligand_resname,
        out_dir,
        args.membrane_core_half_thickness_nm,
    )
    compact_topology, source_numbering_applied = (
        compact_topology_with_source_numbering(
            run_dir,
            prepared["topology"],
            prepared["complex"],
            prepared["protein"],
            prepared["parameterized_bonds"],
            out_dir,
        )
    )
    ligand_xyz, compact_xyz, membrane_center_z = collect_selected_frames(
        trajectory_path,
        selected_frames,
        prepared,
        args.chunk_size,
    )

    cluster_manifest = None
    if not args.skip_clustering:
        _assignments, cluster_summary, cluster_manifest, _labels = (
            cluster_ligand_poses(
                ligand_xyz,
                selected_frames,
                selected_times,
                args.cluster_threshold_A,
                out_dir,
            )
        )
        write_medoid_structures(
            compact_xyz,
            compact_topology,
            cluster_summary,
            out_dir,
            args.max_medoid_structures,
        )

    mmpbsa_rows = []
    mmpbsa_manifest = None
    if prepare_energy:
        if not len(prepared["lipid_heavy"]):
            fail(
                "This endpoint protocol is configured for an implicit membrane, "
                "but no lipid atoms were found."
            )
        topology_result = prepare_amber_topologies(
            run_dir,
            prepared,
            out_dir,
            args.amber_radii,
        )
        trajectory_out, snapshot_frames = write_mmpbsa_trajectory(
            compact_xyz,
            compact_topology,
            membrane_center_z,
            selected_frames,
            selected_times,
            out_dir,
            args.mmpbsa_snapshots,
        )
        input_path = topology_result["directory"] / "mmpbsa.in"
        write_mmpbsa_input(
            input_path,
            len(snapshot_frames),
            args.salt_molar,
            args.solute_dielectric,
            args.membrane_dielectric,
            args.membrane_thickness_A,
        )
        mmpbsa_manifest = {
            "scope": (
                "single-trajectory comparative endpoint MM/PB(GB)SA; "
                "configurational entropy not included"
            ),
            "trajectory": str(trajectory_out),
            "snapshots": int(len(snapshot_frames)),
            "window_ns": [
                float(selected_times[0]),
                float(selected_times[-1]),
            ],
            "models": {
                "GB": {
                    "igb": 5,
                    "salt_molar": args.salt_molar,
                    "membrane": "not represented; sensitivity model only",
                },
                "PB": {
                    "implicit_membrane": True,
                    "solute_dielectric": args.solute_dielectric,
                    "membrane_dielectric": args.membrane_dielectric,
                    "solvent_dielectric": 80.0,
                    "membrane_center_z_A": 0.0,
                    "membrane_thickness_A": args.membrane_thickness_A,
                    "salt_molar": args.salt_molar,
                },
            },
            "comparison_policy": (
                "Use replicate means and standard deviations; compare only "
                "the three SGLT2 ligands generated with the same receptor and "
                "parameterization workflow."
            ),
            "limitations": [
                "The result is an endpoint estimate, not an alchemical free energy.",
                "Configurational entropy is omitted.",
                "The GB calculation is aqueous and does not represent the membrane.",
                "The PB calculation replaces explicit lipids with an implicit membrane slab.",
                "Individual trajectory frames are not independent replicates.",
            ],
        }
        write_json(
            topology_result["directory"] / "mmpbsa_manifest.json",
            mmpbsa_manifest,
        )
        if args.run_mmpbsa:
            mmpbsa_rows = run_mmpbsa(
                topology_result["directory"],
                args.mpi_ranks,
                args.mmpbsa_executable,
            )

    summary = {
        "run_dir": str(run_dir),
        "trajectory_frames": frame_count,
        "trajectory_final_time_ns": float(times_ns[-1]),
        "selected_late_window_ns": [
            float(selected_times[0]),
            float(selected_times[-1]),
        ],
        "selected_frames": int(len(selected_frames)),
        "ligand_resname": prepared["topology"].atom(
            int(prepared["ligand"][0])
        ).residue.name,
        "ligand_heavy_atoms": int(len(prepared["ligand_heavy"])),
        "parameterized_bond_source": prepared["bond_source"],
        "alignment": prepared["alignment_metadata"],
        "source_residue_numbering_applied": source_numbering_applied,
        "clustering": cluster_manifest,
        "mmpbsa_prepared": prepare_energy,
        "mmpbsa_completed": bool(mmpbsa_rows),
        "mmpbsa_results": mmpbsa_rows,
    }
    write_json(out_dir / "endpoint_analysis_summary.json", summary)
    print(f"Endpoint outputs: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
