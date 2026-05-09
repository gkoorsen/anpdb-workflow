"""Check whether Amber-prepared production MD inputs are ready.

Usage:
  python scripts/md_check_amber_inputs.py
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "configs" / "md" / "amber_production"


def resolve(path_text: str) -> Path:
    path = Path(path_text)
    if not path.is_absolute():
        path = ROOT / path
    return path


def check_config(config_path: Path) -> dict[str, object]:
    with config_path.open("rb") as handle:
        cfg = tomllib.load(handle)
    result: dict[str, object] = {
        "config": str(config_path.relative_to(ROOT)),
        "target": cfg["run"]["target"],
        "ready": True,
        "checks": [],
    }
    checks: list[dict[str, object]] = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "ok": ok, "detail": detail})
        if not ok:
            result["ready"] = False

    prmtop = resolve(cfg["input"]["prmtop"])
    inpcrd = resolve(cfg["input"]["inpcrd"])
    add("prmtop_exists", prmtop.exists(), str(prmtop))
    add("inpcrd_exists", inpcrd.exists(), str(inpcrd))
    add("preparation_notes_present", bool(cfg["input"].get("preparation_notes")), "Document cofactor chemistry and prep route")

    result["checks"] = checks
    return result


def main() -> int:
    configs = sorted(CONFIG_DIR.glob("*.toml"))
    results = [check_config(path) for path in configs]
    ready = all(bool(item["ready"]) for item in results)
    print(json.dumps({"ready": ready, "configs": results}, indent=2))
    return 0 if ready else 1


if __name__ == "__main__":
    raise SystemExit(main())
