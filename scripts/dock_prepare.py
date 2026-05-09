"""Prepare receptors + ligands for AutoDock Vina docking on the three lead targets.

Receptors:
  - CYP1B1 (PDB 4I8V), cocrystal ligand BHF (alpha-naphthoflavone), cofactor HEM
  - SGLT2  (PDB 7VSI), cocrystal ligand 7R3 (sotagliflozin)
  - MAO-B  (PDB 2V5Z), cocrystal ligand SAG (safinamide), cofactor FAD

For each receptor:
  1. Extract chain A only (drop dimer / homotetramer copies)
  2. Drop waters and lipids
  3. Extract the cocrystal ligand → use its centre of mass as grid-box centre
  4. Convert the cleaned receptor (with cofactor) to PDBQT via openbabel
  5. Convert each shortlist ligand SMILES to 3D PDBQT via meeko

Output:  output/docking/receptors/{4I8V,7VSI,2V5Z}_receptor.pdbqt
         output/docking/receptors/{4I8V,7VSI,2V5Z}_box.txt
         output/docking/ligands/{compound_id}.pdbqt
         output/docking/dock_jobs.tsv  (target × ligand grid for the runner)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
RECEPTORS = DOCK / "receptors"
LIGANDS = DOCK / "ligands"
LIGANDS.mkdir(parents=True, exist_ok=True)

NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"
CONSENSUS = ROOT / "output" / "chembl_nn" / "consensus_pidgin_chembl_nn_ad60.tsv"

OBABEL = "/opt/homebrew/bin/obabel"
MK_PREP = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/mk_prepare_ligand.py"
PY_AUTODOCK = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/python"

TARGETS = {
    # uniprot, pdb, cocrystal-resname, cofactor-keep, chain
    "CYP1B1": {"uniprot": "Q16678", "pdb": "4I8V", "ligand_resn": "BHF", "keep_het": ["HEM"], "chain": "A"},
    "SGLT2":  {"uniprot": "P31639", "pdb": "7VSI", "ligand_resn": "7R3", "keep_het": [],      "chain": "A"},
    "MAO-B":  {"uniprot": "P27338", "pdb": "2V5Z", "ligand_resn": "SAG", "keep_het": ["FAD"], "chain": "A"},
}

BOX_PADDING = 8.0  # Å around cocrystal ligand centre


def split_pdb(pdb_path: Path, chain: str, ligand_resn: str, keep_het: list[str]):
    """Return (receptor_lines, ligand_lines, ligand_coords) — chain-A only."""
    receptor_lines = []
    ligand_lines = []
    ligand_coords = []
    with open(pdb_path) as fh:
        for line in fh:
            tag = line[:6].strip()
            if tag in {"HEADER", "TITLE", "REMARK", "CRYST1"}:
                receptor_lines.append(line)
                continue
            if tag in {"ATOM", "HETATM"}:
                ch = line[21]
                resn = line[17:20].strip()
                if ch != chain:
                    continue
                if resn == "HOH":
                    continue
                if resn == ligand_resn:
                    ligand_lines.append(line)
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    ligand_coords.append((x, y, z))
                    continue
                if tag == "HETATM":
                    if resn not in keep_het:
                        continue
                receptor_lines.append(line)
            elif tag in {"TER", "END"}:
                receptor_lines.append(line)
    return receptor_lines, ligand_lines, ligand_coords


def write_pdb(lines, path: Path):
    with open(path, "w") as fh:
        fh.writelines(lines)


def receptor_to_pdbqt(in_pdb: Path, out_pdbqt: Path):
    """Use openbabel to add hydrogens at pH 7 and write PDBQT."""
    cmd = [
        OBABEL, str(in_pdb), "-O", str(out_pdbqt),
        "-xh",       # add hydrogens (rigid receptor)
        "-xr",       # treat as rigid (no torsions)
        "-p", "7.4", # protonation pH
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  obabel STDERR: {r.stderr}", file=sys.stderr)


def ligand_to_pdbqt(smi: str, cid: str, out_pdbqt: Path) -> bool:
    """Generate 3D coords + meeko-prepare a ligand from SMILES."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) < 0:
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) < 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass

    sdf_path = out_pdbqt.with_suffix(".sdf")
    Chem.SDWriter(str(sdf_path)).write(mol)

    cmd = [PY_AUTODOCK, MK_PREP, "-i", str(sdf_path), "-o", str(out_pdbqt)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    sdf_path.unlink(missing_ok=True)
    return r.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0


def main() -> int:
    # --- Receptors ---
    boxes: dict[str, dict] = {}
    for tname, t in TARGETS.items():
        pdb = t["pdb"]
        in_pdb = RECEPTORS / f"{pdb}.pdb"
        clean = RECEPTORS / f"{pdb}_clean.pdb"
        ligref = RECEPTORS / f"{pdb}_cocrystal_{t['ligand_resn']}.pdb"
        rec_pdbqt = RECEPTORS / f"{pdb}_receptor.pdbqt"
        print(f"\n=== {tname} ({pdb}) ===")
        rec_lines, lig_lines, lig_coords = split_pdb(
            in_pdb, t["chain"], t["ligand_resn"], t["keep_het"])
        write_pdb(rec_lines, clean)
        write_pdb(lig_lines, ligref)
        if not lig_coords:
            print(f"  WARN: no cocrystal ligand atoms found for {t['ligand_resn']}")
            continue
        coords = np.array(lig_coords)
        cen = coords.mean(axis=0)
        ext = coords.max(axis=0) - coords.min(axis=0) + 2 * BOX_PADDING
        boxes[tname] = {
            "pdb": pdb,
            "center_x": round(cen[0], 3),
            "center_y": round(cen[1], 3),
            "center_z": round(cen[2], 3),
            "size_x":   round(max(ext[0], 22.0), 3),
            "size_y":   round(max(ext[1], 22.0), 3),
            "size_z":   round(max(ext[2], 22.0), 3),
        }
        print(f"  receptor atoms : {sum(1 for l in rec_lines if l.startswith(('ATOM','HETATM')))}")
        print(f"  cocrystal atoms: {len(lig_coords)}")
        print(f"  box centre     : {boxes[tname]['center_x']:.2f} {boxes[tname]['center_y']:.2f} {boxes[tname]['center_z']:.2f}")
        print(f"  box size       : {boxes[tname]['size_x']:.1f} x {boxes[tname]['size_y']:.1f} x {boxes[tname]['size_z']:.1f}")

        receptor_to_pdbqt(clean, rec_pdbqt)
        if rec_pdbqt.exists() and rec_pdbqt.stat().st_size > 0:
            print(f"  receptor PDBQT : {rec_pdbqt.name} ({rec_pdbqt.stat().st_size//1024} KB)")
        else:
            print(f"  receptor PDBQT FAILED")

    pd.DataFrame.from_dict(boxes, orient="index").reset_index(names="target").to_csv(
        RECEPTORS / "boxes.csv", index=False)

    # --- Ligands: 35 unique compounds from consensus set for the three targets ---
    cons = pd.read_csv(CONSENSUS, sep="\t", dtype=str, keep_default_na=False)
    target_uniprots = {t["uniprot"]: tname for tname, t in TARGETS.items()}
    cons["target"] = cons["uniprot"].map(target_uniprots)
    sub = cons[cons["target"].notna()].copy()
    print(f"\n=== Ligands ===")
    print(f"  Consensus pairs across 3 targets: {len(sub)}")
    print(f"  Distinct compounds: {sub['compound_id'].nunique()}")

    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    smi_map = dict(zip(novel["molecule_id"], novel["std_smiles"].fillna("")))

    n_ok = 0
    n_fail = 0
    for cid in sorted(sub["compound_id"].unique()):
        smi = smi_map.get(cid, "")
        if not smi:
            n_fail += 1
            continue
        out_pdbqt = LIGANDS / f"{cid}.pdbqt"
        if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            n_ok += 1
            continue
        if ligand_to_pdbqt(smi, cid, out_pdbqt):
            n_ok += 1
        else:
            n_fail += 1
    print(f"  Prepared PDBQTs: {n_ok}  failed: {n_fail}")

    # --- Job grid ---
    jobs = []
    for _, r in sub.iterrows():
        cid = r["compound_id"]; tname = r["target"]
        if not (LIGANDS / f"{cid}.pdbqt").exists():
            continue
        if tname not in boxes:
            continue
        jobs.append({
            "compound_id": cid,
            "target": tname,
            "uniprot": r["uniprot"],
            "pdb": boxes[tname]["pdb"],
            "ligand_pdbqt":   str(LIGANDS / f"{cid}.pdbqt"),
            "receptor_pdbqt": str(RECEPTORS / f"{boxes[tname]['pdb']}_receptor.pdbqt"),
            **{k: boxes[tname][k] for k in ("center_x","center_y","center_z","size_x","size_y","size_z")},
        })
    pd.DataFrame(jobs).to_csv(DOCK / "dock_jobs.tsv", sep="\t", index=False)
    print(f"\n  Wrote {len(jobs)} dock jobs to dock_jobs.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
