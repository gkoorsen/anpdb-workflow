"""Run the corrected per-replicate and publication MD analysis workflow.

The manuscript set contains three replicates of four corrected OpenFF systems
plus the unaffected MAO-B/Amber system. Use ``--scope corrected`` on the GPU
machine when only the 12 corrected trajectories are present.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Run:
    run_dir: Path
    ligand_resname: str
    corrected: bool


RUNS = [
    *[
        Run(
            ROOT / f"md_runs/production/{system}/rep{replicate}",
            "UNK",
            True,
        )
        for system in (
            "sglt2_mol13144",
            "sglt2_mol13733",
            "sglt2_mol15088",
            "oprk1_mol16614",
        )
        for replicate in (1, 2, 3)
    ],
    Run(
        ROOT / "md_runs/production/maob_mol14056_amber/rep1_rerun2",
        "LIG",
        False,
    ),
    Run(ROOT / "md_runs/production/maob_mol14056_amber/rep2", "LIG", False),
    Run(ROOT / "md_runs/production/maob_mol14056_amber/rep3", "LIG", False),
]


def require_run_inputs(run: Run) -> None:
    required = [
        run.run_dir / "equilibrated.pdb",
        run.run_dir / "production.dcd",
        run.run_dir / "production.log",
        run.run_dir / "config.toml",
        run.run_dir / "run_manifest.json",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        details = "\n".join(f"  - {path}" for path in missing)
        raise SystemExit(f"Missing analysis inputs for {run.run_dir}:\n{details}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scope",
        choices=("all", "corrected"),
        default="all",
        help="Analyze all 15 manuscript runs or only the 12 corrected reruns.",
    )
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--extra-stride", type=int, default=10)
    parser.add_argument("--write-imaged", action="store_true")
    parser.add_argument("--skip-publication", action="store_true")
    args = parser.parse_args()

    selected = [run for run in RUNS if args.scope == "all" or run.corrected]
    for run in selected:
        require_run_inputs(run)

    analyzer = ROOT / "scripts" / "md_analyze_production.py"
    for index, run in enumerate(selected, start=1):
        print(f"[{index}/{len(selected)}] Analyzing {run.run_dir.relative_to(ROOT)}", flush=True)
        command = [
            sys.executable,
            str(analyzer),
            "--run-dir",
            str(run.run_dir),
            "--ligand-resname",
            run.ligand_resname,
            "--stride",
            str(args.stride),
        ]
        if args.write_imaged:
            command.append("--write-imaged")
        subprocess.run(command, cwd=ROOT, check=True)

    if not args.skip_publication:
        if args.scope != "all":
            raise SystemExit(
                "Publication aggregation requires --scope all. "
                "Use --skip-publication when analyzing only corrected runs."
            )
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "md_publication_analysis.py"),
                "--extra-stride",
                str(args.extra_stride),
            ],
            cwd=ROOT,
            check=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
