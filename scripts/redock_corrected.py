"""Redo the redock control with cocrystal SMILES extracted directly from
each PDB cocrystal extract (not from PubChem look-ups).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
RECEPTORS = DOCK / "receptors"
VAL = DOCK / "validation"

VINA = "/usr/local/bin/vina"
PY_AUTODOCK = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/python"
MK_PREP = "/Users/gerritkoorsen/opt/anaconda3/envs/autodock/bin/mk_prepare_ligand.py"
OBABEL = "/opt/homebrew/bin/obabel"

# SMILES extracted directly from the cocrystal PDB by OpenBabel
COCRYSTAL = {
    "CYP1B1": ("4I8V", "BHF", "alpha-naphthoflavone (BHF)",
               "c1c(=O)c2ccc3ccccc3c2oc1c1ccccc1"),
    "SGLT2":  ("7VSI", "7R3", "empagliflozin (7R3)",
               "[C@H]1([C@H]([C@@H]([C@@H](CO)O[C@H]1c1ccc(c(Cc2ccc(cc2)O[C@@H]2COCC2)c1)Cl)O)O)O"),
    "MAO-B":  ("2V5Z", "SAG", "safinamide Schiff base (SAG)",
               "Fc1cc(ccc1)COc1ccc(cc1)/C=N/[C@@H](C)C(=O)N"),
}

AUTODOCK_TO_ELEMENT = {
    "A": "C", "C": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "HD": "H", "H": "H",
    "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I",
    "P": "P",
}


def normalise_element(raw: str, atom_name: str = "") -> str:
    raw = raw.strip()
    if raw in AUTODOCK_TO_ELEMENT:
        return AUTODOCK_TO_ELEMENT[raw]
    name = atom_name.strip()
    if name:
        letters = ""
        for ch in name:
            if ch.isalpha():
                letters += ch
            else:
                break
        if len(letters) >= 2 and letters[:2].title() in AUTODOCK_TO_ELEMENT:
            return AUTODOCK_TO_ELEMENT[letters[:2].title()]
        if letters[:1].upper() in AUTODOCK_TO_ELEMENT:
            return AUTODOCK_TO_ELEMENT[letters[:1].upper()]
    cleaned = "".join(c for c in raw if c.isalpha()).capitalize()
    return cleaned


def parse_pdb_heavy(path: Path):
    atoms = []
    in_model = True
    seen_model = False
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            seen_model = True
            in_model = (line.split()[1].strip() == "1")
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if not seen_model:
            in_model = True
        if not in_model:
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            raw_elem = line[76:78] if len(line) > 78 else ""
            atom_name = line[12:16] if len(line) > 16 else ""
            elem = normalise_element(raw_elem, atom_name)
            if elem == "H" or not elem:
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append((elem, x, y, z))
            except ValueError:
                continue
    return atoms


def greedy_rmsd(crystal, posed):
    used = [False] * len(posed)
    sq = []
    matched = 0
    for ce, cx, cy, cz in crystal:
        best = None; best_idx = -1
        for i, (pe, px, py, pz) in enumerate(posed):
            if used[i] or pe != ce:
                continue
            d = (cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2
            if best is None or d < best:
                best = d; best_idx = i
        if best is None:
            continue  # skip unmatched (asymmetric truncation in Schiff base etc.)
        used[best_idx] = True
        sq.append(best)
        matched += 1
    if matched < 3:
        return None, matched
    return float(np.sqrt(np.mean(sq))), matched


def ligand_to_pdbqt(smi: str, out_pdbqt: Path) -> bool:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return False
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) < 0:
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) < 0: return False
    try: AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception: pass
    sdf = out_pdbqt.with_suffix(".sdf")
    Chem.SDWriter(str(sdf)).write(mol)
    r = subprocess.run([PY_AUTODOCK, MK_PREP, "-i", str(sdf), "-o", str(out_pdbqt)],
                       capture_output=True, text=True)
    sdf.unlink(missing_ok=True)
    return out_pdbqt.exists() and out_pdbqt.stat().st_size > 0


def vina_dock(receptor: Path, ligand: Path, box: dict, out_pdbqt: Path):
    cmd = [VINA, "--receptor", str(receptor), "--ligand", str(ligand),
           "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
           "--center_z", str(box["center_z"]),
           "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]),
           "--size_z", str(box["size_z"]),
           "--out", str(out_pdbqt),
           "--exhaustiveness", "16", "--num_modes", "9", "--seed", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    aff = None
    if out_pdbqt.exists():
        for line in out_pdbqt.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                try: aff = float(line.split()[3]); break
                except: pass
    return aff


def main() -> int:
    boxes = pd.read_csv(RECEPTORS / "boxes.csv").set_index("target").to_dict("index")

    rows = []
    for tname, (pdb, resn, name, smi) in COCRYSTAL.items():
        print(f"=== {tname} ({pdb} / {resn}) ===")
        print(f"   SMILES: {smi}")
        rec = RECEPTORS / f"{pdb}_receptor.pdbqt"
        cocrys_pdb = RECEPTORS / f"{pdb}_cocrystal_{resn}.pdb"
        lig_pdbqt = VAL / f"redock2_{tname}_{resn}.pdbqt"
        out_pdbqt = VAL / f"redock2_{tname}_{resn}_out.pdbqt"
        if not ligand_to_pdbqt(smi, lig_pdbqt):
            print(f"   FAILED ligand prep")
            continue
        aff = vina_dock(rec, lig_pdbqt, boxes[tname], out_pdbqt)
        print(f"   affinity: {aff} kcal/mol")

        c = parse_pdb_heavy(cocrys_pdb)
        p = parse_pdb_heavy(out_pdbqt)
        rmsd, matched = greedy_rmsd(c, p)
        rows.append({
            "target": tname, "pdb": pdb, "ligand": name, "resname": resn,
            "best_affinity": aff,
            "n_xtal_heavy": len(c), "n_pose_heavy": len(p),
            "n_matched": matched,
            "rmsd_absolute_A": round(rmsd, 2) if rmsd is not None else None,
        })
        if rmsd is not None:
            print(f"   xtal={len(c)} pose={len(p)} matched={matched} RMSD(abs)={rmsd:.2f} A")
        else:
            print(f"   xtal={len(c)} pose={len(p)} matched={matched} RMSD=N/A")

    df = pd.DataFrame(rows)
    df.to_csv(VAL / "redock_corrected.tsv", sep="\t", index=False)
    print("\nFinal:")
    print(df.to_string(index=False))
    n_pass = sum(1 for _, r in df.iterrows()
                 if r["rmsd_absolute_A"] is not None and r["rmsd_absolute_A"] < 2.0)
    n_marginal = sum(1 for _, r in df.iterrows()
                     if r["rmsd_absolute_A"] is not None and 2.0 <= r["rmsd_absolute_A"] < 3.0)
    print(f"\n{n_pass}/{len(df)} cocrystal redocks pass <2 A RMSD; {n_marginal} marginal (2-3 A)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
