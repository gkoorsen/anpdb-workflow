"""Config-driven OpenMM production MD runner.

This script is intentionally stricter than scripts/md_run.py. It is meant for
long GPU production runs and records enough metadata to support reproducible
methods.

Typical usage:
  python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml
  python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --resume
  python scripts/md_production.py --config configs/md/dev/apo_water_10ps.toml --allow-non-production
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform as py_platform
import shutil
import subprocess
import sys
import time
import tomllib
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PRODUCTION_GATES = {
    "CYP1B1": {
        "required_heterogens": {"HEM"},
        "expected_cofactor_residues": {"HEM"},
        "requires_membrane": False,
        "note": "CYP1B1 production runs should retain and parameterise heme.",
    },
    "MAO-B": {
        "required_heterogens": {"FAD"},
        "expected_cofactor_residues": {"FAD"},
        "requires_membrane": False,
        "note": "MAO-B production runs should retain and parameterise FAD.",
    },
    "SGLT2": {
        "required_heterogens": set(),
        "expected_cofactor_residues": set(),
        "requires_membrane": True,
        "note": "SGLT2 production runs should use a membrane-oriented structure and lipid bilayer.",
    },
}


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


def read_resnames(pdb_path: Path) -> set[str]:
    resnames: set[str] = set()
    with pdb_path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                resnames.add(line[17:20].strip())
    return resnames


def ffxml_residue_names(path: Path) -> set[str]:
    tree = ET.parse(path)
    names: set[str] = set()
    for elem in tree.getroot().iter():
        if elem.tag.rsplit("}", 1)[-1] == "Residue" and "name" in elem.attrib:
            names.add(elem.attrib["name"])
    return names


def validate_production_gates(config: dict[str, Any], receptor_pdb: Path | None, allow_non_production: bool) -> None:
    target = str(require(config, "run.target"))
    require_production = bool(deep_get(config, "production.require_production_ready", True))
    if not require_production:
        if not allow_non_production:
            fail("Config disables production gates. Re-run with --allow-non-production for dev/apo-water tests.")
        return

    gate = PRODUCTION_GATES.get(target)
    if gate is None:
        fail(f"No production gate defined for target: {target}")

    errors: list[str] = []
    if gate["requires_membrane"]:
        environment = deep_get(config, "system.environment", "solvent")
        membrane_enabled = bool(deep_get(config, "system.membrane.enabled", False))
        oriented = bool(deep_get(config, "input.receptor_is_membrane_oriented", False))
        if environment != "membrane" or not membrane_enabled:
            errors.append("SGLT2 requires system.environment: membrane and system.membrane.enabled: true.")
        if not oriented:
            errors.append("SGLT2 requires input.receptor_is_membrane_oriented: true after OPM/CHARMM-GUI/manual orientation.")

    required_heterogens = set(gate["required_heterogens"])
    if required_heterogens:
        cofactor_policy = deep_get(config, "system.cofactors.policy", "strip")
        parameter_files = deep_get(config, "system.cofactors.parameter_files", [])
        if cofactor_policy != "require":
            errors.append(f"{target} requires system.cofactors.policy: require.")
        if not parameter_files:
            errors.append(f"{target} requires cofactor parameter files for {', '.join(sorted(required_heterogens))}.")
        if receptor_pdb and receptor_pdb.exists():
            present = read_resnames(receptor_pdb)
            missing = required_heterogens - present
            if missing:
                errors.append(f"Receptor is missing required heterogens: {', '.join(sorted(missing))}.")

    if errors:
        lines = [f"Production gates failed for {target}: {gate['note']}"]
        lines.extend(f"- {err}" for err in errors)
        fail("\n".join(lines))


def select_inputs(config: dict[str, Any], root: Path) -> dict[str, Any]:
    compound_id = deep_get(config, "input.compound_id")
    target = str(require(config, "run.target"))
    pdb = deep_get(config, "input.pdb")

    dock_results = resolve_path(root, deep_get(config, "input.dock_results", "output/docking/dock_results.tsv"))
    if not compound_id:
        if dock_results is None or not dock_results.exists():
            fail("input.compound_id is not set and dock_results file is missing.")
        best: dict[str, str] | None = None
        best_affinity: float | None = None
        with dock_results.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("target") != target or not row.get("best_affinity"):
                    continue
                affinity = float(row["best_affinity"])
                if best_affinity is None or affinity < best_affinity:
                    best = row
                    best_affinity = affinity
        if best is None:
            fail(f"No docking results found for target {target}")
        compound_id = str(best["compound_id"])
        pdb = str(best["pdb"])

    if not pdb:
        fail("input.pdb must be set when input.compound_id is set.")

    receptor_pdb = resolve_path(root, deep_get(config, "input.receptor_pdb"))
    if receptor_pdb is None:
        receptor_pdb = root / "output" / "docking" / "receptors" / f"{pdb}_clean.pdb"

    ligand_pdbqt = resolve_path(root, deep_get(config, "input.ligand_pdbqt"))
    if ligand_pdbqt is None:
        ligand_pdbqt = root / "output" / "docking" / "results" / f"{compound_id}_{pdb}_out.pdbqt"

    smiles_table = resolve_path(root, deep_get(config, "input.ligand_smiles_table", "output/anpdb_truly_novel_std.csv"))
    return {
        "compound_id": compound_id,
        "target": target,
        "pdb": pdb,
        "receptor_pdb": receptor_pdb,
        "ligand_pdbqt": ligand_pdbqt,
        "smiles_table": smiles_table,
    }


def check_input_paths(config: dict[str, Any], root: Path, paths: dict[str, Any]) -> None:
    missing: list[str] = []
    for key, path in paths.items():
        if not (key.endswith("_pdb") or key.endswith("_pdbqt") or key.endswith("_table")):
            continue
        if isinstance(path, Path) and not path.exists():
            missing.append(str(path))

    target = str(require(config, "run.target"))
    gate = PRODUCTION_GATES.get(target, {})
    expected_cofactor_residues = set(gate.get("expected_cofactor_residues", set()))
    cofactor_residues: set[str] = set()

    for path_value in deep_get(config, "system.cofactors.parameter_files", []) or []:
        path = resolve_path(root, path_value)
        if path is not None and not path.exists():
            missing.append(str(path))
        elif path is not None:
            try:
                cofactor_residues.update(ffxml_residue_names(path))
            except Exception as exc:
                missing.append(f"{path} is not a parseable OpenMM ffxml file: {exc}")

    if expected_cofactor_residues:
        missing_residues = expected_cofactor_residues - cofactor_residues
        if missing_residues:
            missing.append(
                "cofactor XML files do not define expected residue templates: "
                + ", ".join(sorted(missing_residues))
            )

    if missing:
        fail("Input readiness checks failed:\n" + "\n".join(f"- {path}" for path in missing))


def collect_input_issues(config: dict[str, Any], root: Path, paths: dict[str, Any]) -> list[str]:
    try:
        check_input_paths(config, root, paths)
    except SystemExit as exc:
        return [line.removeprefix("- ") for line in str(exc).splitlines()[1:]]
    return []


def extract_pdbqt_model(pdbqt_path: Path, mode: int, out_pdbqt: Path) -> None:
    keep: list[str] = []
    in_mode = False
    saw = 0
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL"):
            saw += 1
            in_mode = saw == mode
            if in_mode:
                keep.append(line)
            continue
        if line.startswith("ENDMDL"):
            if in_mode:
                keep.append(line)
                break
            in_mode = False
            continue
        if in_mode:
            keep.append(line)
    if not keep:
        fail(f"Could not extract mode {mode} from {pdbqt_path}")
    out_pdbqt.write_text("\n".join(keep) + "\n")


def convert_pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path, obabel: str) -> None:
    result = subprocess.run(
        [obabel, str(pdbqt_path), "-O", str(pdb_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not pdb_path.exists() or pdb_path.stat().st_size == 0:
        fail(
            f"Open Babel conversion failed ({obabel}):\n"
            f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
        )


def load_ligand_smiles(smiles_table: Path, compound_id: str, id_column: str, smiles_column: str) -> str:
    with smiles_table.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or id_column not in reader.fieldnames or smiles_column not in reader.fieldnames:
            fail(f"SMILES table must contain columns {id_column!r} and {smiles_column!r}: {smiles_table}")
        for row in reader:
            if row.get(id_column) != compound_id:
                continue
            smiles = str(row.get(smiles_column, ""))
            if not smiles:
                fail(f"Empty SMILES for {compound_id}")
            return smiles.replace("[O][Na]", "[O-]").replace("[O][K]", "[O-]")
    fail(f"No SMILES found for {compound_id} in {smiles_table}")


def prepare_openff_ligand(config: dict[str, Any], pose_pdb: Path, smiles: str):
    from md_run import pose_to_openff_mol
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    mode = deep_get(config, "ligand.pose_mode", "graft")
    charge_model = deep_get(config, "ligand.charge_model", "openff-gnn-am1bcc-0.1.0-rc.3.pt")
    lig_mol = pose_to_openff_mol(pose_pdb, smiles, mode=mode)
    lig_mol.assign_partial_charges(charge_model, toolkit_registry=NAGLToolkitWrapper())
    return lig_mol


def write_receptor_for_policy(config: dict[str, Any], receptor_pdb: Path, out_pdb: Path) -> None:
    policy = deep_get(config, "system.cofactors.policy", "strip")
    strip_resnames = set(deep_get(config, "system.cofactors.strip_resnames", ["HEM", "HEME", "FAD"]))
    if policy == "strip":
        with receptor_pdb.open() as fin, out_pdb.open("w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() in strip_resnames:
                    continue
                fout.write(line)
    elif policy == "require":
        shutil.copyfile(receptor_pdb, out_pdb)
    else:
        fail(f"Unsupported system.cofactors.policy: {policy}")


def ps_to_steps(ps: float, timestep_fs: float) -> int:
    return int(round(ps * 1000.0 / timestep_fs))


def ns_to_steps(ns: float, timestep_fs: float) -> int:
    return ps_to_steps(ns * 1000.0, timestep_fs)


def interval_steps(ps: float, timestep_fs: float) -> int:
    return max(1, ps_to_steps(ps, timestep_fs))


def build_simulation(config: dict[str, Any], root: Path, paths: dict[str, Any], run_dir: Path):
    from openmm import MonteCarloBarostat, Platform, unit
    from openmm.app import ForceField, HBonds, Modeller, PDBFile, PME, Simulation
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from pdbfixer import PDBFixer
    from openmm import LangevinMiddleIntegrator

    mode = int(deep_get(config, "input.docking_mode", deep_get(config, "input.pose_mode", 1)))
    obabel = deep_get(config, "tools.obabel", "obabel")
    pose_pdbqt = run_dir / f"{paths['compound_id']}_mode{mode}.pdbqt"
    pose_pdb = run_dir / f"{paths['compound_id']}_mode{mode}.pdb"
    extract_pdbqt_model(paths["ligand_pdbqt"], mode, pose_pdbqt)
    convert_pdbqt_to_pdb(pose_pdbqt, pose_pdb, obabel)

    smiles = load_ligand_smiles(
        paths["smiles_table"],
        paths["compound_id"],
        deep_get(config, "input.ligand_id_column", "molecule_id"),
        deep_get(config, "input.ligand_smiles_column", "std_smiles"),
    )
    lig_mol = prepare_openff_ligand(config, pose_pdb, smiles)

    receptor_policy_pdb = run_dir / "receptor_policy_input.pdb"
    receptor_fixed_pdb = run_dir / "receptor_fixed.pdb"
    write_receptor_for_policy(config, paths["receptor_pdb"], receptor_policy_pdb)

    fixer = PDBFixer(filename=str(receptor_policy_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if deep_get(config, "system.cofactors.policy", "strip") == "strip":
        fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens(float(deep_get(config, "system.ph", 7.4)))
    PDBFile.writeFile(fixer.topology, fixer.positions, receptor_fixed_pdb.open("w"), keepIds=True)

    rec = PDBFile(str(receptor_fixed_pdb))
    modeller = Modeller(rec.topology, rec.positions)
    modeller.add(lig_mol.to_topology().to_openmm(), lig_mol.conformers[0].to_openmm())

    ff_files = list(deep_get(config, "forcefield.files", ["amber/ff14SB.xml", "amber/tip3p_standard.xml"]))
    cofactor_files = deep_get(config, "system.cofactors.parameter_files", []) or []
    ff = ForceField(*(ff_files + [str(resolve_path(root, p)) for p in cofactor_files]))
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[lig_mol])
    ff.registerTemplateGenerator(smirnoff.generator)

    environment = deep_get(config, "system.environment", "solvent")
    ionic_strength = float(deep_get(config, "system.ionic_strength_molar", 0.15)) * unit.molar
    if environment == "membrane":
        membrane = deep_get(config, "system.membrane", {})
        modeller.addMembrane(
            ff,
            lipidType=membrane.get("lipid_type", "POPC"),
            membraneCenterZ=float(membrane.get("center_z_nm", 0.0)) * unit.nanometer,
            minimumPadding=float(membrane.get("minimum_padding_nm", 1.0)) * unit.nanometer,
            ionicStrength=ionic_strength,
            positiveIon=deep_get(config, "system.positive_ion", "Na+"),
            negativeIon=deep_get(config, "system.negative_ion", "Cl-"),
        )
    elif environment == "solvent":
        modeller.addSolvent(
            ff,
            model=deep_get(config, "system.water_model", "tip3p"),
            padding=float(deep_get(config, "system.padding_nm", 1.0)) * unit.nanometer,
            ionicStrength=ionic_strength,
            positiveIon=deep_get(config, "system.positive_ion", "Na+"),
            negativeIon=deep_get(config, "system.negative_ion", "Cl-"),
        )
    else:
        fail(f"Unsupported system.environment: {environment}")

    solvated_pdb = run_dir / "system_solvated.pdb"
    PDBFile.writeFile(modeller.topology, modeller.positions, solvated_pdb.open("w"))

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=float(deep_get(config, "simulation.nonbonded_cutoff_nm", 1.0)) * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=bool(deep_get(config, "simulation.remove_cm_motion", False)),
        hydrogenMass=float(deep_get(config, "simulation.hydrogen_mass_amu", 4.0)) * unit.amu,
    )
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

    platform_name = deep_get(config, "platform.name", "CUDA")
    platform = Platform.getPlatformByName(platform_name)
    properties = {}
    if platform_name == "CUDA":
        properties = {
            "DeviceIndex": str(deep_get(config, "platform.device_index", "0")),
            "Precision": str(deep_get(config, "platform.precision", "mixed")),
        }
    sim = Simulation(modeller.topology, system, integrator, platform, properties)
    sim.context.setPositions(modeller.positions)
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


def run_protocol(config: dict[str, Any], root: Path, paths: dict[str, Any], run_dir: Path, resume: bool) -> None:
    import openmm
    from openmm import unit
    from openmm.app import PDBFile

    sim, platform, properties = build_simulation(config, root, paths, run_dir)
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
        append = True
    else:
        min_iters = int(deep_get(config, "simulation.minimization_max_iterations", 5000))
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
        append = False

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
            "target": paths["target"],
            "compound_id": paths["compound_id"],
            "pdb": paths["pdb"],
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
    parser.add_argument("--allow-non-production", action="store_true")
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = load_config(config_path)
    root = repo_root_from_config(config, config_path)
    paths = select_inputs(config, root)
    validate_production_gates(config, paths.get("receptor_pdb"), args.allow_non_production)

    run_dir = resolve_path(root, require(config, "run.output_dir"))
    assert run_dir is not None
    plan = {
        "config": str(config_path),
        "root": str(root),
        "run_dir": str(run_dir),
        "paths": {key: str(value) if isinstance(value, Path) else value for key, value in paths.items()},
        "production_ns": deep_get(config, "simulation.production_ns"),
        "platform": deep_get(config, "platform.name", "CUDA"),
    }
    if args.dry_run:
        input_issues = collect_input_issues(config, root, paths)
        plan["input_ready"] = not input_issues
        plan["input_issues"] = input_issues
        print(json.dumps(plan, indent=2))
        return 0

    check_input_paths(config, root, paths)
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_path, run_dir / f"config{config_path.suffix}")
    write_json(run_dir / "plan.json", plan)

    started = time.time()
    write_json(
        run_dir / "run_started.json",
        {
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "pid": os.getpid(),
            "argv": sys.argv,
        },
    )
    run_protocol(config, root, paths, run_dir, resume=args.resume)
    print(f"Finished in {(time.time() - started) / 3600:.2f} h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
