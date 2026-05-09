"""Check whether the ignored production MD input bundle is ready.

Usage:
  python scripts/md_check_inputs.py
"""

from __future__ import annotations

import json
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "md" / "production"
EXPECTED_COFACTOR_RESIDUES = {
    "CYP1B1": {"HEM"},
    "MAO-B": {"FAD"},
}


def resnames_in_pdb(path: Path) -> set[str]:
    resnames: set[str] = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                resnames.add(line[17:20].strip())
    return resnames


def resolve(path_text: str) -> Path:
    path = Path(path_text)
    if not path.is_absolute():
        path = ROOT / path
    return path


def xml_is_parseable(path: Path) -> bool:
    try:
        ET.parse(path)
    except Exception:
        return False
    return True


def ffxml_residue_names(path: Path) -> set[str]:
    tree = ET.parse(path)
    names: set[str] = set()
    for elem in tree.getroot().iter():
        if elem.tag.rsplit("}", 1)[-1] == "Residue" and "name" in elem.attrib:
            names.add(elem.attrib["name"])
    return names


def check_config(config_path: Path) -> dict[str, object]:
    with config_path.open("rb") as handle:
        cfg = tomllib.load(handle)
    target = cfg["run"]["target"]
    result: dict[str, object] = {
        "config": str(config_path.relative_to(ROOT)),
        "target": target,
        "ready": True,
        "checks": [],
    }
    checks: list[dict[str, object]] = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "ok": ok, "detail": detail})
        if not ok:
            result["ready"] = False

    receptor = resolve(cfg["input"]["receptor_pdb"])
    add("receptor_exists", receptor.exists(), str(receptor))
    if receptor.exists():
        resnames = resnames_in_pdb(receptor)
        if target == "CYP1B1":
            add("heme_present", "HEM" in resnames, "HEM must be present in CYP1B1 receptor")
        if target == "MAO-B":
            add("fad_present", "FAD" in resnames, "FAD must be present in MAO-B receptor")

    ligand_pose = resolve(cfg["input"]["ligand_pdbqt"])
    add("ligand_pose_exists", ligand_pose.exists(), str(ligand_pose))

    smiles = resolve(cfg["input"]["ligand_smiles_table"])
    add("smiles_table_exists", smiles.exists(), str(smiles))

    if target == "SGLT2":
        oriented = bool(cfg["input"].get("receptor_is_membrane_oriented", False))
        add("membrane_orientation_declared", oriented, "Config must declare receptor_is_membrane_oriented = true")
        add("membrane_environment", cfg["system"].get("environment") == "membrane", "system.environment must be membrane")
        add("membrane_enabled", bool(cfg["system"].get("membrane", {}).get("enabled", False)), "system.membrane.enabled must be true")

    cofactor_residues: set[str] = set()
    for cofactor_path in cfg["system"].get("cofactors", {}).get("parameter_files", []):
        xml_path = resolve(cofactor_path)
        add("cofactor_xml_exists", xml_path.exists(), str(xml_path))
        if xml_path.exists():
            parseable = xml_is_parseable(xml_path)
            add("cofactor_xml_parseable", parseable, str(xml_path))
            if parseable:
                cofactor_residues.update(ffxml_residue_names(xml_path))

    expected_residues = EXPECTED_COFACTOR_RESIDUES.get(target, set())
    if expected_residues:
        missing = expected_residues - cofactor_residues
        add(
            "cofactor_templates_match_target",
            not missing,
            "Expected OpenMM ffxml Residue templates: " + ", ".join(sorted(expected_residues)),
        )

    result["checks"] = checks
    return result


def main() -> int:
    configs = sorted(CONFIG_DIR.glob("*.toml"))
    results = [check_config(path) for path in configs]
    ready = all(bool(item["ready"]) for item in results)
    payload = {"ready": ready, "configs": results}
    print(json.dumps(payload, indent=2))
    return 0 if ready else 1


if __name__ == "__main__":
    raise SystemExit(main())
