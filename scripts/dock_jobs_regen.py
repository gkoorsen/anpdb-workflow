"""Regenerate dock_jobs.tsv after the 5 ligand failures were recovered."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
LIGANDS = DOCK / "ligands"
RECEPTORS = DOCK / "receptors"
CONSENSUS = ROOT / "output" / "chembl_nn" / "consensus_pidgin_chembl_nn_ad60.tsv"

TARGETS = {
    "Q16678": ("CYP1B1", "4I8V"),
    "P31639": ("SGLT2", "7VSI"),
    "P27338": ("MAO-B", "2V5Z"),
}


def main() -> int:
    boxes = pd.read_csv(RECEPTORS / "boxes.csv").set_index("target").to_dict("index")
    cons = pd.read_csv(CONSENSUS, sep="\t", dtype=str, keep_default_na=False)
    cons = cons[cons["uniprot"].isin(TARGETS)].copy()
    cons["target"] = cons["uniprot"].map(lambda u: TARGETS[u][0])
    cons["pdb"] = cons["uniprot"].map(lambda u: TARGETS[u][1])

    jobs = []
    for _, r in cons.iterrows():
        cid = r["compound_id"]; tname = r["target"]
        ligfile = LIGANDS / f"{cid}.pdbqt"
        if not ligfile.exists() or ligfile.stat().st_size == 0:
            continue
        b = boxes[tname]
        jobs.append({
            "compound_id": cid,
            "target": tname,
            "uniprot": r["uniprot"],
            "pdb": r["pdb"],
            "ligand_pdbqt":   str(ligfile),
            "receptor_pdbqt": str(RECEPTORS / f"{r['pdb']}_receptor.pdbqt"),
            "center_x": b["center_x"], "center_y": b["center_y"], "center_z": b["center_z"],
            "size_x":   b["size_x"],   "size_y":   b["size_y"],   "size_z":   b["size_z"],
        })
    df = pd.DataFrame(jobs)
    df.to_csv(DOCK / "dock_jobs.tsv", sep="\t", index=False)
    print(f"  Wrote {len(df)} dock jobs")
    print(df.groupby("target").size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
