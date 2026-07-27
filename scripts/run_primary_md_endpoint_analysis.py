"""Run late-pose clustering and endpoint energies for the corrected MD set.

Clustering is performed for every available corrected SGLT2 and OPRK1
trajectory.  MM/PB(GB)SA preparation/execution is restricted to SGLT2 so that
energetic comparisons remain within one receptor and parameterization workflow.
Lightweight outputs are copied into the tracked manuscript results directory.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = (
    ROOT
    / "results"
    / "md_analysis"
    / "manuscript_v2"
    / "endpoint_analysis"
)
SYSTEMS = (
    "sglt2_mol13144",
    "sglt2_mol13733",
    "sglt2_mol15088",
    "oprk1_mol16614",
)
LIGHTWEIGHT_FILES = (
    "alignment_core_atoms.csv",
    "cluster_assignments.csv",
    "cluster_summary.csv",
    "clustering_manifest.json",
    "dominant_cluster_medoid.pdb",
    "endpoint_analysis_summary.json",
    "residue_number_map.csv",
)
LIGHTWEIGHT_MMPBSA_FILES = (
    "FINAL_RESULTS_MMPBSA.csv",
    "FINAL_RESULTS_MMPBSA.dat",
    "mmpbsa.in",
    "mmpbsa_command.json",
    "mmpbsa_manifest.json",
    "mmpbsa_summary.csv",
    "mmpbsa_summary.json",
    "snapshot_frames.csv",
    "topology_conversion_manifest.json",
)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def required_inputs(
    run_dir: Path,
    require_solvated_topology: bool,
) -> list[Path]:
    paths = [
        run_dir / "equilibrated.pdb",
        run_dir / "production.dcd",
        run_dir / "production.log",
        run_dir / "config.toml",
        run_dir / "system.xml",
    ]
    if require_solvated_topology:
        paths.append(run_dir / "system_solvated.pdb")
    return paths


def run_directories() -> list[tuple[str, int, Path]]:
    return [
        (
            system,
            replicate,
            ROOT / "md_runs" / "production" / system / f"rep{replicate}",
        )
        for system in SYSTEMS
        for replicate in (1, 2, 3)
    ]


def copy_if_present(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True


def export_lightweight_outputs(
    system: str,
    replicate: int,
    run_dir: Path,
) -> None:
    source = run_dir / "endpoint_analysis"
    destination = RESULTS_ROOT / system / f"rep{replicate}"
    destination.mkdir(parents=True, exist_ok=True)
    for name in LIGHTWEIGHT_FILES:
        copy_if_present(source / name, destination / name)
    mmpbsa_source = source / "mmpbsa"
    mmpbsa_destination = destination / "mmpbsa"
    for name in LIGHTWEIGHT_MMPBSA_FILES:
        copy_if_present(
            mmpbsa_source / name,
            mmpbsa_destination / name,
        )


def aggregate_results(completed: list[tuple[str, int, Path]]) -> None:
    cluster_rows = []
    energy_rows = []
    for system, replicate, _run_dir in completed:
        result_dir = RESULTS_ROOT / system / f"rep{replicate}"
        endpoint_summary_path = result_dir / "endpoint_analysis_summary.json"
        if endpoint_summary_path.exists():
            endpoint = json.loads(endpoint_summary_path.read_text())
            clustering = endpoint.get("clustering")
            if clustering:
                cluster_rows.append(
                    {
                        "system": system,
                        "replicate": replicate,
                        "selected_frames": clustering["selected_frames"],
                        "clusters": clustering["clusters"],
                        "cluster_threshold_A": clustering[
                            "distance_threshold_A"
                        ],
                        "dominant_cluster_fraction": clustering[
                            "dominant_cluster_fraction"
                        ],
                        "dominant_medoid_frame_index": clustering[
                            "dominant_cluster_medoid_frame_index"
                        ],
                        "dominant_medoid_time_ns": clustering[
                            "dominant_cluster_medoid_time_ns"
                        ],
                        "mean_silhouette": clustering.get(
                            "mean_silhouette"
                        ),
                    }
                )
        energy_summary_path = result_dir / "mmpbsa" / "mmpbsa_summary.json"
        if energy_summary_path.exists():
            energy = json.loads(energy_summary_path.read_text())
            for row in energy.get("models", []):
                energy_rows.append(
                    {
                        "system": system,
                        "replicate": replicate,
                        **row,
                    }
                )

    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    if cluster_rows:
        pd.DataFrame(cluster_rows).sort_values(
            ["system", "replicate"]
        ).to_csv(
            RESULTS_ROOT / "clustering_replicate_summary.csv",
            index=False,
        )
    if energy_rows:
        energy = pd.DataFrame(energy_rows).sort_values(
            ["system", "replicate", "model"]
        )
        energy.to_csv(
            RESULTS_ROOT / "mmpbsa_replicate_results.csv",
            index=False,
        )
        grouped_rows = []
        for (system, model), group in energy.groupby(
            ["system", "model"],
            sort=True,
        ):
            replicate_values = group["delta_endpoint_kcal_mol"]
            grouped_rows.append(
                {
                    "system": system,
                    "model": model,
                    "replicates": int(len(group)),
                    "replicate_mean_kcal_mol": float(
                        replicate_values.mean()
                    ),
                    "replicate_sd_kcal_mol": (
                        float(replicate_values.std(ddof=1))
                        if len(group) > 1
                        else None
                    ),
                    "replicate_values_kcal_mol": ";".join(
                        f"{value:.6f}"
                        for value in replicate_values
                    ),
                }
            )
        pd.DataFrame(grouped_rows).to_csv(
            RESULTS_ROOT / "mmpbsa_complex_summary.csv",
            index=False,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scope",
        choices=("completed", "all"),
        default="completed",
        help=(
            "Skip incomplete runs or require all 12 corrected trajectories."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("cluster", "prepare", "run"),
        default="run",
        help=(
            "Cluster only; cluster and prepare MMPBSA.py inputs; or execute "
            "the complete workflow."
        ),
    )
    parser.add_argument("--late-window-ns", type=float, default=20.0)
    parser.add_argument("--cluster-stride", type=int, default=1)
    parser.add_argument("--cluster-threshold-A", type=float, default=2.0)
    parser.add_argument("--chunk-size", type=int, default=25)
    parser.add_argument("--mmpbsa-snapshots", type=int, default=100)
    parser.add_argument("--mpi-ranks", type=int, default=1)
    parser.add_argument("--mmpbsa-executable")
    parser.add_argument(
        "--skip-clustering",
        action="store_true",
        help="Reuse existing cluster outputs while preparing/running energies.",
    )
    args = parser.parse_args()

    require_energy_inputs = args.mode in {"prepare", "run"}
    selected = []
    missing_rows = []
    for system, replicate, run_dir in run_directories():
        need_system = require_energy_inputs and system.startswith("sglt2_")
        missing = [
            path
            for path in required_inputs(
                run_dir,
                require_solvated_topology=need_system,
            )
            if not path.exists()
        ]
        if missing:
            missing_rows.append(
                {
                    "system": system,
                    "replicate": replicate,
                    "run_dir": str(run_dir),
                    "missing": [str(path) for path in missing],
                }
            )
            if args.scope == "all":
                raise SystemExit(
                    f"Run is incomplete: {run_dir}\n"
                    + "\n".join(f"  - {path}" for path in missing)
                )
            print(f"Skipping incomplete run: {run_dir.relative_to(ROOT)}")
            continue
        selected.append((system, replicate, run_dir))

    if not selected:
        raise SystemExit("No completed corrected trajectories were found.")

    analyzer = ROOT / "scripts" / "md_endpoint_analysis.py"
    completed = []
    for index, (system, replicate, run_dir) in enumerate(
        selected,
        start=1,
    ):
        print(
            f"[{index}/{len(selected)}] {system}/rep{replicate}",
            flush=True,
        )
        command = [
            sys.executable,
            str(analyzer),
            "--run-dir",
            str(run_dir),
            "--ligand-resname",
            "UNK",
            "--late-window-ns",
            str(args.late_window_ns),
            "--cluster-stride",
            str(args.cluster_stride),
            "--cluster-threshold-A",
            str(args.cluster_threshold_A),
            "--chunk-size",
            str(args.chunk_size),
            "--mmpbsa-snapshots",
            str(args.mmpbsa_snapshots),
            "--mpi-ranks",
            str(args.mpi_ranks),
        ]
        if args.skip_clustering:
            command.append("--skip-clustering")
        if system.startswith("sglt2_") and args.mode in {"prepare", "run"}:
            command.append("--prepare-mmpbsa")
        if system.startswith("sglt2_") and args.mode == "run":
            command.append("--run-mmpbsa")
        if args.mmpbsa_executable:
            command.extend(
                ["--mmpbsa-executable", args.mmpbsa_executable]
            )
        subprocess.run(command, cwd=ROOT, check=True)
        export_lightweight_outputs(system, replicate, run_dir)
        completed.append((system, replicate, run_dir))

    aggregate_results(completed)
    write_json(
        RESULTS_ROOT / "batch_manifest.json",
        {
            "scope": args.scope,
            "mode": args.mode,
            "completed_runs": [
                {
                    "system": system,
                    "replicate": replicate,
                    "run_dir": str(run_dir),
                }
                for system, replicate, run_dir in completed
            ],
            "skipped_runs": missing_rows,
            "mmpbsa_scope": (
                "SGLT2 only; within-target comparison across Mol_13144, "
                "Mol_13733 and Mol_15088"
            ),
            "replication_policy": (
                "Independent replicate means are the unit of replication; "
                "individual snapshots are not treated as independent samples."
            ),
        },
    )
    print(f"Tracked lightweight results: {RESULTS_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
