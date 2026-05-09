"""Run production MD from curated Amber topology/coordinate files.

Use this runner for targets whose cofactors are best handled during system
preparation, for example CYP1B1 heme and covalent MAO-B FAD. The input Amber
files must already contain the protein, ligand, cofactors, solvent/ions, and
any covalent cofactor patches.

Typical usage:
  python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml
  python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --equilibrate-only
  python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --resume
"""

from __future__ import annotations

import argparse
import json
import os
import platform as py_platform
import shutil
import sys
import time
import tomllib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def fail(message: str) -> None:
    raise SystemExit(message)


def load_config(path: Path) -> dict[str, Any]:
    if path.suffix.lower() == ".toml":
        with path.open("rb") as handle:
            data = tomllib.load(handle)
    else:
        try:
            import yaml
        except Exception as exc:
            fail(f"PyYAML is required to read {path.name}. Use a .toml config or install pyyaml. ({exc})")
        with path.open() as handle:
            data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        fail(f"Config is empty or invalid: {path}")
    return data


def deep_get(data: dict[str, Any], dotted: str, default: Any = None) -> Any:
    cur: Any = data
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def require(data: dict[str, Any], dotted: str) -> Any:
    value = deep_get(data, dotted)
    if value in (None, ""):
        fail(f"Missing required config key: {dotted}")
    return value


def repo_root_from_config(config: dict[str, Any], config_path: Path) -> Path:
    root_value = deep_get(config, "project.root", ".")
    root = Path(root_value)
    if not root.is_absolute():
        root = (config_path.parent / root).resolve()
        if not (root / "scripts").exists():
            root = Path.cwd().resolve()
    return root


def resolve_path(root: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = root / path
    return path


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def ps_to_steps(ps: float, timestep_fs: float) -> int:
    return int(round(ps * 1000.0 / timestep_fs))


def ns_to_steps(ns: float, timestep_fs: float) -> int:
    return ps_to_steps(ns * 1000.0, timestep_fs)


def interval_steps(ps: float, timestep_fs: float) -> int:
    return max(1, ps_to_steps(ps, timestep_fs))


def select_inputs(config: dict[str, Any], root: Path) -> dict[str, Any]:
    prmtop = resolve_path(root, str(require(config, "input.prmtop")))
    inpcrd = resolve_path(root, str(require(config, "input.inpcrd")))
    assert prmtop is not None and inpcrd is not None
    return {
        "compound_id": deep_get(config, "input.compound_id"),
        "target": str(require(config, "run.target")),
        "pdb": deep_get(config, "input.pdb"),
        "prmtop": prmtop,
        "inpcrd": inpcrd,
        "prepared_by": deep_get(config, "input.prepared_by"),
        "preparation_notes": deep_get(config, "input.preparation_notes"),
    }


def collect_input_issues(paths: dict[str, Any]) -> list[str]:
    issues: list[str] = []
    for key in ("prmtop", "inpcrd"):
        path = paths[key]
        if isinstance(path, Path) and not path.exists():
            issues.append(str(path))
    return issues


def platform_properties(config: dict[str, Any]) -> tuple[str, dict[str, str]]:
    platform_name = str(deep_get(config, "platform.name", "CUDA"))
    properties: dict[str, str] = {}
    if platform_name == "CUDA":
        properties = {
            "DeviceIndex": str(deep_get(config, "platform.device_index", "0")),
            "Precision": str(deep_get(config, "platform.precision", "mixed")),
        }
    return platform_name, properties


def build_simulation(config: dict[str, Any], paths: dict[str, Any]):
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME, Simulation

    inpcrd = AmberInpcrdFile(str(paths["inpcrd"]))
    prmtop = AmberPrmtopFile(str(paths["prmtop"]), periodicBoxVectors=inpcrd.boxVectors)
    hydrogen_mass_amu = float(deep_get(config, "simulation.hydrogen_mass_amu", 4.0))
    hydrogen_mass = hydrogen_mass_amu * unit.amu if hydrogen_mass_amu > 0 else None
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=float(deep_get(config, "simulation.nonbonded_cutoff_nm", 1.0)) * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=bool(deep_get(config, "simulation.remove_cm_motion", False)),
        hydrogenMass=hydrogen_mass,
    )
    if bool(deep_get(config, "simulation.add_barostat", True)):
        system.addForce(
            MonteCarloBarostat(
                float(deep_get(config, "simulation.pressure_atm", 1.0)) * unit.atmospheres,
                float(deep_get(config, "simulation.temperature_K", 310.0)) * unit.kelvin,
                int(deep_get(config, "simulation.barostat_interval_steps", 25)),
            )
        )

    integrator = LangevinMiddleIntegrator(
        float(deep_get(config, "simulation.temperature_K", 310.0)) * unit.kelvin,
        float(deep_get(config, "simulation.friction_per_ps", 1.0)) / unit.picosecond,
        float(deep_get(config, "simulation.timestep_fs", 2.0)) * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(int(deep_get(config, "run.seed", 1)))

    platform_name, properties = platform_properties(config)
    platform = Platform.getPlatformByName(platform_name)
    sim = Simulation(prmtop.topology, system, integrator, platform, properties)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    return sim, platform, properties


def add_reporters(config: dict[str, Any], sim, run_dir: Path, append: bool, production_final_step: int) -> None:
    from openmm.app import CheckpointReporter, DCDReporter, StateDataReporter

    timestep_fs = float(deep_get(config, "simulation.timestep_fs", 2.0))
    sim.reporters.append(
        DCDReporter(
            str(run_dir / "production.dcd"),
            interval_steps(float(deep_get(config, "output.trajectory_interval_ps", 50.0)), timestep_fs),
            append=append,
        )
    )
    sim.reporters.append(
        StateDataReporter(
            str(run_dir / "production.log"),
            interval_steps(float(deep_get(config, "output.state_interval_ps", 50.0)), timestep_fs),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            speed=True,
            remainingTime=True,
            totalSteps=production_final_step,
            append=append,
        )
    )
    sim.reporters.append(
        CheckpointReporter(
            str(run_dir / "production.chk"),
            interval_steps(float(deep_get(config, "output.checkpoint_interval_ps", 500.0)), timestep_fs),
        )
    )


def run_protocol(config: dict[str, Any], paths: dict[str, Any], run_dir: Path, resume: bool, equilibrate_only: bool) -> None:
    import openmm
    from openmm import unit
    from openmm.app import PDBFile

    sim, platform, properties = build_simulation(config, paths)
    timestep_fs = float(deep_get(config, "simulation.timestep_fs", 2.0))
    production_steps = ns_to_steps(float(deep_get(config, "simulation.production_ns", 100.0)), timestep_fs)

    checkpoint = run_dir / "production.chk"
    production_plan_path = run_dir / "production_plan.json"
    resumed = resume and checkpoint.exists()
    if resumed:
        if not production_plan_path.exists():
            fail(f"Cannot resume without {production_plan_path}")
        production_plan = json.loads(production_plan_path.read_text())
        production_start_step = int(production_plan["production_start_step"])
        production_final_step = int(production_plan["production_final_step"])
        sim.loadCheckpoint(str(checkpoint))
        append = (run_dir / "production.dcd").exists() and (run_dir / "production.log").exists()
    else:
        min_iters = int(deep_get(config, "simulation.minimization_max_iterations", 5000))
        if min_iters > 0:
            sim.minimizeEnergy(maxIterations=min_iters)
        state = sim.context.getState(getPositions=True)
        PDBFile.writeFile(sim.topology, state.getPositions(), (run_dir / "minimized.pdb").open("w"))

        sim.context.setVelocitiesToTemperature(
            float(deep_get(config, "simulation.temperature_K", 310.0)) * unit.kelvin,
            int(deep_get(config, "run.seed", 1)),
        )
        nvt_ps = float(deep_get(config, "simulation.equilibration.nvt_ps", 250.0))
        npt_ps = float(deep_get(config, "simulation.equilibration.npt_ps", 1000.0))
        if nvt_ps > 0:
            sim.step(ps_to_steps(nvt_ps, timestep_fs))
        if npt_ps > 0:
            sim.step(ps_to_steps(npt_ps, timestep_fs))
        state = sim.context.getState(getPositions=True, getVelocities=True)
        PDBFile.writeFile(sim.topology, state.getPositions(), (run_dir / "equilibrated.pdb").open("w"))
        sim.saveState(str(run_dir / "equilibrated_state.xml"))
        production_start_step = int(sim.currentStep)
        production_final_step = production_start_step + production_steps
        write_json(
            production_plan_path,
            {
                "production_start_step": production_start_step,
                "production_steps": production_steps,
                "production_final_step": production_final_step,
                "production_ns": deep_get(config, "simulation.production_ns"),
                "timestep_fs": timestep_fs,
            },
        )
        sim.saveCheckpoint(str(checkpoint))
        write_json(
            run_dir / "equilibration_manifest.json",
            {
                "completed_utc": datetime.now(timezone.utc).isoformat(),
                "python": sys.version.split()[0],
                "openmm": openmm.version.version,
                "host_platform": py_platform.platform(),
                "input_system_type": "amber",
                "target": paths["target"],
                "compound_id": paths["compound_id"],
                "pdb": paths["pdb"],
                "prmtop": str(paths["prmtop"]),
                "inpcrd": str(paths["inpcrd"]),
                "prepared_by": paths.get("prepared_by"),
                "preparation_notes": paths.get("preparation_notes"),
                "run_name": deep_get(config, "run.name"),
                "replicate": deep_get(config, "run.replicate"),
                "seed": deep_get(config, "run.seed"),
                "platform": platform.getName(),
                "platform_properties": properties,
                "minimization_max_iterations": min_iters,
                "nvt_ps": nvt_ps,
                "npt_ps": npt_ps,
                "timestep_fs": timestep_fs,
                "production_start_step": production_start_step,
                "production_checkpoint": str(checkpoint),
                "next_command": "rerun this config with --resume to start production from production.chk",
            },
        )
        append = False

    if equilibrate_only:
        print(f"Equilibration complete. Resume production with --resume from {checkpoint}")
        return

    add_reporters(config, sim, run_dir, append=append, production_final_step=production_final_step)

    remaining = production_final_step - int(sim.currentStep)
    if remaining <= 0:
        print("Production already complete according to checkpoint.", file=sys.stderr)
    else:
        sim.step(remaining)

    state = sim.context.getState(getPositions=True, getVelocities=True)
    PDBFile.writeFile(sim.topology, state.getPositions(), (run_dir / "final.pdb").open("w"))
    sim.saveState(str(run_dir / "final_state.xml"))

    write_json(
        run_dir / "run_manifest.json",
        {
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.split()[0],
            "openmm": openmm.version.version,
            "host_platform": py_platform.platform(),
            "input_system_type": "amber",
            "target": paths["target"],
            "compound_id": paths["compound_id"],
            "pdb": paths["pdb"],
            "prmtop": str(paths["prmtop"]),
            "inpcrd": str(paths["inpcrd"]),
            "prepared_by": paths.get("prepared_by"),
            "preparation_notes": paths.get("preparation_notes"),
            "run_name": deep_get(config, "run.name"),
            "replicate": deep_get(config, "run.replicate"),
            "seed": deep_get(config, "run.seed"),
            "platform": platform.getName(),
            "platform_properties": properties,
            "production_ns": deep_get(config, "simulation.production_ns"),
            "timestep_fs": timestep_fs,
            "production_start_step": production_start_step,
            "production_final_step": production_final_step,
            "resumed": resumed,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--equilibrate-only", action="store_true", help="Run minimization/equilibration, write a production-start checkpoint, then stop.")
    args = parser.parse_args()
    if args.equilibrate_only and args.resume:
        fail("Use --equilibrate-only without --resume. After inspection, run production with --resume.")

    config_path = args.config.resolve()
    config = load_config(config_path)
    root = repo_root_from_config(config, config_path)
    paths = select_inputs(config, root)
    run_dir = resolve_path(root, require(config, "run.output_dir"))
    assert run_dir is not None

    plan = {
        "config": str(config_path),
        "root": str(root),
        "run_dir": str(run_dir),
        "input_system_type": "amber",
        "paths": {key: str(value) if isinstance(value, Path) else value for key, value in paths.items()},
        "production_ns": deep_get(config, "simulation.production_ns"),
        "equilibrate_only": args.equilibrate_only,
        "platform": deep_get(config, "platform.name", "CUDA"),
    }
    input_issues = collect_input_issues(paths)
    plan["input_ready"] = not input_issues
    plan["input_issues"] = input_issues
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    if input_issues:
        fail("Input readiness checks failed:\n" + "\n".join(f"- {issue}" for issue in input_issues))

    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_path, run_dir / f"config{config_path.suffix}")
    write_json(run_dir / "plan.json", plan)
    write_json(
        run_dir / "run_started.json",
        {
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "pid": os.getpid(),
            "argv": sys.argv,
        },
    )
    started = time.time()
    run_protocol(config, paths, run_dir, resume=args.resume, equilibrate_only=args.equilibrate_only)
    print(f"Finished in {(time.time() - started) / 3600:.2f} h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
