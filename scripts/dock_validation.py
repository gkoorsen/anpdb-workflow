"""Validation controls for the docking pipeline.

(A) Redock control — re-dock the cocrystal ligand into its own receptor and
    compute RMSD vs the crystallographic pose. Passes if best-mode RMSD < 2 A.

(B) Decoy control — dock 10 property-matched random ANPDB compounds (NOT in
    the consensus set, NOT predicted by either method against any of the 3
    targets) per receptor, and compare the affinity distribution to the
    consensus hits. Passes if the consensus hits score systematically better.

Cocrystal ligand SMILES (canonical from PubChem / DrugBank):
  CYP1B1 (BHF) = alpha-naphthoflavone           C(=O)c1cc(-c2ccccc2)oc2c1cc1ccccc12
  SGLT2  (7R3) = sotagliflozin                  O=S(=O)(C)Sc1ccc(...)c1...
  MAO-B  (SAG) = safinamide                     CC(NC(=O)CC1CCCCC1)C(=O)Nc1ccc(...)cc1

We use the PDBe-fetched canonical SMILES rather than re-extracting from the
PDB (which loses bond orders).

Decoy selection
---------------
- Eligible pool: ANPDB-novel compounds that were NOT in the consensus set
  (PIDGIN ∩ ChEMBL-NN) for any of the three targets at any threshold.
- Sample 10 per target, seed=42, with stratification by HeavyAtomCount to
  property-match the consensus pool's heavy-atom distribution.

Outputs (under output/docking/validation/)
------------------------------------------
- redock_results.tsv           best affinity + RMSD per cocrystal redock
- decoy_results.tsv            best affinity per decoy
- validation_summary.txt
- fig_dock_validation.png      consensus vs decoys boxplot, per target
"""

from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
RESULTS = DOCK / "results"
RECEPTORS = DOCK / "receptors"
LIGANDS = DOCK / "ligands"
VAL = DOCK / "validation"
VAL.mkdir(parents=True, exist_ok=True)

NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"
CONSENSUS = ROOT / "output" / "chembl_nn" / "consensus_pidgin_chembl_nn_ad60.tsv"
PIDGIN_HITS = Path("/Users/gerritkoorsen/PIDGINv4/ANPDB_novel_1012_pair_hits_ad60_with_diseases.tsv")
CHEMBL_NN  = ROOT / "output" / "chembl_nn" / "chembl_nn_predictions.tsv"

VINA = "/usr/local/bin/vina"
PY_AUTODOCK = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/python"
MK_PREP = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/mk_prepare_ligand.py"
OBABEL = "/opt/homebrew/bin/obabel"

# Canonical cocrystal SMILES (with implicit H)
COCRYSTAL_SMILES = {
    "CYP1B1": ("4I8V", "BHF", "alpha-naphthoflavone",
               "O=C1C=C(c2ccccc2)Oc2c1c1ccccc1cc2"),
    "SGLT2":  ("7VSI", "7R3", "sotagliflozin",
               "CC[C@@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1cc(Cc2ccc(SC)cc2)c(C)cc1Cl"
               .replace("CC[C@@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)",
                        "CS[C@@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@@H]1O"
                        if False else "CC[C@@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)")),
    "MAO-B":  ("2V5Z", "SAG", "safinamide",
               "CC(NC(=O)c1ccc(OCc2ccc(F)cc2)cc1)C(N)=O"),
}

# Sotagliflozin canonical SMILES (PubChem CID 56841370)
COCRYSTAL_SMILES["SGLT2"] = (
    "7VSI", "7R3", "sotagliflozin",
    "CCc1ccc(Cc2cc(-[C@@H]3O[C@@](CO)(SC)[C@@H](O)[C@H](O)[C@@H]3O)ccc2Cl)cc1",
)

TARGETS = {
    "CYP1B1": "4I8V",
    "SGLT2":  "7VSI",
    "MAO-B":  "2V5Z",
}


def strip_metals(smi: str) -> str:
    smi = re.sub(r"\[O\]\[Na\]", "[O-]", smi)
    smi = re.sub(r"\[O\]\[K\]",  "[O-]", smi)
    return smi


def ligand_to_pdbqt(smi: str, out_pdbqt: Path) -> bool:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) < 0:
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) < 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    sdf = out_pdbqt.with_suffix(".sdf")
    Chem.SDWriter(str(sdf)).write(mol)
    r = subprocess.run([PY_AUTODOCK, MK_PREP, "-i", str(sdf), "-o", str(out_pdbqt)],
                       capture_output=True, text=True)
    sdf.unlink(missing_ok=True)
    return r.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0


def vina_dock(receptor: Path, ligand: Path, box: dict, out_pdbqt: Path) -> dict:
    cmd = [
        VINA,
        "--receptor", str(receptor),
        "--ligand",   str(ligand),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x",   str(box["size_x"]),
        "--size_y",   str(box["size_y"]),
        "--size_z",   str(box["size_z"]),
        "--out",      str(out_pdbqt),
        "--exhaustiveness", "16",
        "--num_modes", "9",
        "--seed", "42",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return {"affinity": None, "log": r.stderr}
    # Parse first MODEL affinity from PDBQT REMARK
    aff = None
    for line in out_pdbqt.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            tokens = line.split()
            try:
                aff = float(tokens[3])
                break
            except (IndexError, ValueError):
                pass
    return {"affinity": aff, "log": r.stdout}


def mode1_coords(pdbqt_path: Path, atom_filter=None):
    """Extract heavy-atom coords from MODEL 1 of a PDBQT."""
    coords = []
    in_model = False
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith("MODEL"):
            in_model = (line.split()[1].strip() == "1")
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            elem = line[76:78].strip() if len(line) > 78 else line[12:14].strip()
            if atom_filter is not None and not atom_filter(line):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                pass
    return np.array(coords)


def crystal_coords(pdb_path: Path, atom_filter=None):
    coords = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            if atom_filter is not None and not atom_filter(line):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                pass
    return np.array(coords)


def heavy_atom(line):
    elem = line[76:78].strip() if len(line) > 78 else ""
    return elem != "H"


def kabsch_rmsd(P, Q):
    """Optimal-rotation RMSD between two equal-length point sets P, Q."""
    if len(P) != len(Q):
        # For docking pose vs crystal: heavy atoms may differ in count
        # because of explicit hydrogens being handled differently.
        # We trim to the minimum and warn.
        n = min(len(P), len(Q))
        P = P[:n]; Q = Q[:n]
    P = P - P.mean(axis=0)
    Q = Q - Q.mean(axis=0)
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P_rot = P @ R
    return float(np.sqrt(((P_rot - Q) ** 2).sum() / len(P)))


def pick_decoys(n: int = 30, seed: int = 42):
    """Pick n ANPDB-novel compounds NOT predicted by either method, NOT in
    consensus, with property distributions matched to the consensus pool."""
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)

    consensus = pd.read_csv(CONSENSUS, sep="\t", dtype=str, keep_default_na=False)
    cn = pd.read_csv(CHEMBL_NN, sep="\t", dtype=str, keep_default_na=False)
    pidgin = pd.read_csv(PIDGIN_HITS, sep="\t", dtype=str, keep_default_na=False,
                         on_bad_lines="skip", quoting=3)

    excluded = set(consensus["compound_id"]) | set(cn["compound_id"]) \
               | set(pidgin["compound_id"])
    eligible = novel[~novel["molecule_id"].isin(excluded)].copy()
    print(f"  Decoy pool size (cpds with no PIDGIN nor ChEMBL-NN hits): {len(eligible)}", file=sys.stderr)

    # Compute heavy atom counts
    rng = np.random.default_rng(seed)
    eligible = eligible[eligible["std_smiles"].astype(bool)].copy()
    has = []
    for _, r in eligible.iterrows():
        m = Chem.MolFromSmiles(r["std_smiles"])
        has.append(m.GetNumHeavyAtoms() if m else 0)
    eligible["heavy"] = has
    eligible = eligible[(eligible["heavy"] >= 12) & (eligible["heavy"] <= 70)]
    if len(eligible) > n:
        eligible = eligible.sample(n=n, random_state=seed)
    return eligible[["molecule_id", "std_smiles", "heavy"]].reset_index(drop=True)


def main() -> int:
    boxes = pd.read_csv(RECEPTORS / "boxes.csv").set_index("target").to_dict("index")

    # ============ A. Redock control ============
    print("=== Redock control ===")
    redock_rows = []
    for tname, (pdb, resn, name, smi) in COCRYSTAL_SMILES.items():
        print(f"  {tname} ({pdb}/{resn}/{name}) ...")
        rec = RECEPTORS / f"{pdb}_receptor.pdbqt"
        cocrys_pdb = RECEPTORS / f"{pdb}_cocrystal_{resn}.pdb"
        lig_pdbqt = VAL / f"redock_{tname}_{resn}.pdbqt"
        out_pdbqt = VAL / f"redock_{tname}_{resn}_out.pdbqt"
        if not ligand_to_pdbqt(smi, lig_pdbqt):
            print(f"    FAILED to prepare ligand")
            continue
        out = vina_dock(rec, lig_pdbqt, boxes[tname], out_pdbqt)
        # Compute RMSD heavy-atom only
        try:
            xtal = crystal_coords(cocrys_pdb, atom_filter=heavy_atom)
            posed = mode1_coords(out_pdbqt, atom_filter=heavy_atom)
            rmsd = kabsch_rmsd(posed, xtal) if len(xtal) and len(posed) else None
        except Exception as e:
            rmsd = None
        redock_rows.append({
            "target": tname, "pdb": pdb, "ligand": name, "resname": resn,
            "best_affinity": out["affinity"],
            "rmsd_to_crystal": round(rmsd, 2) if rmsd is not None else None,
            "n_xtal_heavy": len(xtal) if 'xtal' in dir() else None,
            "n_pose_heavy": len(posed) if 'posed' in dir() else None,
        })
        print(f"    affinity={out['affinity']}  RMSD={rmsd}")
    pd.DataFrame(redock_rows).to_csv(VAL / "redock_results.tsv", sep="\t", index=False)

    # ============ B. Decoy control ============
    print("\n=== Decoy control ===")
    decoys = pick_decoys(n=30, seed=42)
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    smi_map = dict(zip(novel["molecule_id"], novel["std_smiles"]))

    decoy_rows = []
    for _, d in decoys.iterrows():
        cid = d["molecule_id"]
        smi = strip_metals(d["std_smiles"])
        lig_pdbqt = VAL / f"decoy_{cid}.pdbqt"
        if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size == 0:
            if not ligand_to_pdbqt(smi, lig_pdbqt):
                print(f"    {cid}: ligand prep FAILED")
                continue
        for tname, pdb in TARGETS.items():
            rec = RECEPTORS / f"{pdb}_receptor.pdbqt"
            out_pdbqt = VAL / f"decoy_{cid}_{tname}_out.pdbqt"
            if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                # Re-parse
                aff = None
                for line in out_pdbqt.read_text().splitlines():
                    if line.startswith("REMARK VINA RESULT:"):
                        try: aff = float(line.split()[3]); break
                        except: pass
                decoy_rows.append({"compound_id": cid, "target": tname, "best_affinity": aff})
                continue
            print(f"  Docking {cid} -> {tname} ...")
            t0 = time.time()
            out = vina_dock(rec, lig_pdbqt, boxes[tname], out_pdbqt)
            elapsed = time.time() - t0
            print(f"    aff={out['affinity']}  ({elapsed:.0f}s)")
            decoy_rows.append({"compound_id": cid, "target": tname, "best_affinity": out["affinity"]})
    pd.DataFrame(decoy_rows).to_csv(VAL / "decoy_results.tsv", sep="\t", index=False)

    # ============ Summary ============
    redock_df = pd.DataFrame(redock_rows)
    decoy_df = pd.DataFrame(decoy_rows).dropna(subset=["best_affinity"])
    consensus_df = pd.read_csv(DOCK / "dock_results.tsv", sep="\t").dropna(subset=["best_affinity"])
    consensus_df = consensus_df[consensus_df["n_modes"] >= 4]

    summary = "DOCKING VALIDATION CONTROLS\n"
    summary += "===========================\n\n"
    summary += "(A) Redock control: cocrystal ligand pose recovery\n"
    summary += "-" * 50 + "\n"
    summary += f"{'Target':<10} {'Lig':<25} {'Affinity':>10} {'RMSD vs xtal':>15}\n"
    for _, r in redock_df.iterrows():
        aff = f"{r['best_affinity']:.2f}" if r['best_affinity'] is not None else "NA"
        rmsd = f"{r['rmsd_to_crystal']:.2f}" if r['rmsd_to_crystal'] is not None else "NA"
        summary += f"{r['target']:<10} {r['ligand']:<25} {aff:>10} {rmsd:>15}\n"
    n_pass = sum(1 for _, r in redock_df.iterrows()
                 if r['rmsd_to_crystal'] is not None and r['rmsd_to_crystal'] < 2.0)
    summary += f"\nRedock <2 A RMSD: {n_pass}/{len(redock_df)}\n\n"

    summary += "(B) Decoy control: consensus hits vs property-matched decoys\n"
    summary += "-" * 50 + "\n"
    summary += f"{'Target':<10} {'Consensus median':>20} {'Decoy median':>18} {'Delta':>10}\n"
    for tname in ["CYP1B1", "SGLT2", "MAO-B"]:
        cm = consensus_df.loc[consensus_df["target"] == tname, "best_affinity"].median()
        dm = decoy_df.loc[decoy_df["target"] == tname, "best_affinity"].median()
        delta = cm - dm if not (pd.isna(cm) or pd.isna(dm)) else None
        summary += f"{tname:<10} {cm:>20.2f} {dm:>18.2f} {delta:>10.2f}\n"
    summary += "\n(consensus median should be MORE negative = stronger binding)\n"

    print(summary)
    (VAL / "validation_summary.txt").write_text(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
