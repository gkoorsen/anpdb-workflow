"""Prepare the ignored local input bundle for production MD.

This script stages the files that can be prepared from this repository's local
docking outputs. It deliberately does not fabricate the hard scientific inputs:

- SGLT2 must still be supplied as an OPM/PPM/CHARMM-GUI-oriented receptor.
- OPRK1 must be supplied in the EncoMPASS/OPM membrane frame with the ligand pose transformed into the same coordinates.
- CYP1B1 heme and MAO-B FAD force-field XML files must still be supplied.

Usage:
  python scripts/md_prepare_inputs.py
  python scripts/md_prepare_inputs.py --force
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
DEST = ROOT / "data" / "md_inputs"


FILES = {
    "smiles": (OUT / "anpdb_truly_novel_std.csv", DEST / "anpdb_truly_novel_std.csv"),
    "cyp1b1_pose": (
        OUT / "docking" / "results" / "Mol_11315_4I8V_out.pdbqt",
        DEST / "poses" / "Mol_11315_4I8V_out.pdbqt",
    ),
    "sglt2_pose": (
        OUT / "docking" / "results" / "Mol_13144_7VSI_out.pdbqt",
        DEST / "poses" / "Mol_13144_7VSI_out.pdbqt",
    ),
    "maob_pose": (
        OUT / "docking" / "results" / "Mol_14056_2V5Z_out.pdbqt",
        DEST / "poses" / "Mol_14056_2V5Z_out.pdbqt",
    ),
    "oprk1_pose": (
        OUT / "docking" / "disease_link_outstanding" / "results" / "Mol_16614_4DJH_out.pdbqt",
        DEST / "poses" / "Mol_16614_4DJH_out.pdbqt",
    ),
    "cyp1b1_receptor": (
        OUT / "docking" / "receptors" / "4I8V_clean.pdb",
        DEST / "receptors" / "4I8V_chainA_heme_prepared.pdb",
    ),
    "maob_receptor": (
        OUT / "docking" / "receptors" / "2V5Z_clean.pdb",
        DEST / "receptors" / "2V5Z_chainA_fad_prepared.pdb",
    ),
    "sglt2_unoriented_reference": (
        OUT / "docking" / "receptors" / "7VSI_clean.pdb",
        DEST / "receptors" / "7VSI_clean_unoriented_reference.pdb",
    ),
    "oprk1_receptor": (
        OUT / "docking" / "disease_link_outstanding" / "receptors" / "4DJH_OPRK1_clean.pdb",
        DEST / "receptors" / "4DJH_OPRK1_clean_unoriented_reference.pdb",
    ),
}


def resnames_in_pdb(path: Path) -> set[str]:
    resnames: set[str] = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                resnames.add(line[17:20].strip())
    return resnames


def copy_file(src: Path, dst: Path, force: bool) -> dict[str, object]:
    if not src.exists():
        return {"source": str(src), "dest": str(dst), "status": "missing_source"}
    if dst.exists() and not force:
        return {"source": str(src), "dest": str(dst), "status": "exists"}
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return {"source": str(src), "dest": str(dst), "status": "copied"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    for subdir in ("poses", "receptors", "cofactors", "amber_systems"):
        (DEST / subdir).mkdir(parents=True, exist_ok=True)

    manifest: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "files": {},
        "production_blockers": [
            "Run scripts/md_fetch_orient_sglt2.py to create the OPM-oriented 7VSI receptor and transformed ligand pose.",
            "Run scripts/md_fetch_orient_oprk1.py to create the EncoMPASS/OPM-oriented 4DJH receptor and transformed Mol_16614 pose.",
            "Inspect OPRK1 equilibrated.pdb carefully because 4DJH is a GPCR crystal construct with a T4 lysozyme fusion segment.",
            "Provide data/md_inputs/cofactors/heme.xml matching HEM in 4I8V_chainA_heme_prepared.pdb.",
            "Provide data/md_inputs/cofactors/fad.xml matching FAD in 2V5Z_chainA_fad_prepared.pdb.",
        ],
    }

    file_status: dict[str, object] = {}
    for label, (src, dst) in FILES.items():
        status = copy_file(src, dst, args.force)
        if status["status"] in {"copied", "exists"} and str(dst).endswith(".pdb"):
            status["resnames"] = sorted(resnames_in_pdb(dst))
        file_status[label] = status
    manifest["files"] = file_status

    readme = DEST / "README.md"
    readme.write_text(
        "# MD Input Bundle\n\n"
        "This directory is intentionally ignored by Git.\n\n"
        "Staged from local outputs:\n\n"
        "- `anpdb_truly_novel_std.csv`\n"
        "- ligand PDBQT docking poses for Mol_11315, Mol_13144, and Mol_14056\n"
        "- ligand PDBQT docking pose for OPRK1 / Mol_16614\n"
        "- CYP1B1 receptor with HEM retained\n"
        "- MAO-B receptor with FAD retained\n"
        "- SGLT2 unoriented reference receptor only\n\n"
        "- OPRK1 4DJH unoriented reference receptor\n\n"
        "Still required before all production runs:\n\n"
        "- `receptors/7VSI_opm_oriented_clean.pdb`\n"
        "- `poses/Mol_13144_7VSI_opm_oriented_out.pdbqt`\n"
        "- `receptors/4DJH_OPRK1_opm_oriented_clean.pdb`\n"
        "- `poses/Mol_16614_4DJH_opm_oriented_out.pdbqt`\n"
        "- inspection of OPRK1 `equilibrated.pdb` after membrane insertion\n"
        "- `cofactors/heme.xml`\n"
        "- `cofactors/fad.xml`\n"
        "- or curated Amber systems under `amber_systems/` for CYP1B1 and MAO-B\n\n"
        "Create the SGLT2-oriented receptor and ligand pose with:\n\n"
        "```bash\n"
        "python scripts/md_fetch_orient_sglt2.py\n"
        "python scripts/md_fetch_orient_oprk1.py\n"
        "```\n",
    )
    (DEST / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
