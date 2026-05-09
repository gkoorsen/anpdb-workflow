"""Pharmacophore-constrained redock analysis for CYP1B1.

Approach:
  1. Find heme Fe atom coordinates in 4I8V receptor.
  2. Find the C4-carbonyl oxygen atom in the crystal-pose alpha-naphthoflavone
     (the chromen-4-one C=O — the canonical CYP-pharmacophore interaction).
  3. Measure crystal Fe ... O(C4=O) distance — the reference constraint.
  4. For each of the 9 docked modes, identify the C4-carbonyl O and measure
     its Fe distance.
  5. Filter modes whose Fe-O distance is within crystal +/- 1.5 A.
  6. Among constraint-satisfying modes, pick the best-scoring.
  7. Report the constrained-best mode's RMSD vs crystal.

This emulates pharmacophore-constrained docking (Korhonen et al. 2007) without
needing smina/Vina-Carb — the constraint is applied as a post-filter.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
VAL = ROOT / "output" / "docking" / "validation"
RECEPTORS = ROOT / "output" / "docking" / "receptors"
RECEPTOR = RECEPTORS / "4I8V_receptor.pdbqt"
COCRYS = RECEPTORS / "4I8V_cocrystal_BHF.pdb"
OBABEL = "/opt/homebrew/bin/obabel"
OBRMS = "/opt/homebrew/bin/obrms"

CONSTRAINT_TOLERANCE_A = 1.5  # +/- A around crystal Fe-O distance


def find_heme_iron(pdbqt: Path):
    """Locate FE in the receptor PDBQT — match by atom name 'FE' (residue name
    may be HEM, UNK, or other depending on OpenBabel's relabelling)."""
    for line in pdbqt.read_text().splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            elem_col = line[76:78].strip() if len(line) > 78 else ""
            if atom_name == "FE" or elem_col.lower() == "fe":
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                return np.array([x, y, z])
    return None


def find_chromen4one_carbonyl_O(pdb_or_pdbqt: Path):
    """Identify the chromen-4-one carbonyl O — atom name 'O1' or similar in BHF.
    For PDBQT we infer from atom name + bond context."""
    # In the 4I8V cocrystal entry, BHF atoms are named like O1, C1...C16
    # The C4 carbonyl O is bonded to C2 (the chromone carbonyl carbon) and is the only sp2 O.
    # In RCSB PDB chemical component ID BHF, atom O1 is the chromenone C=O
    # (https://www.rcsb.org/ligand/BHF). Use that.
    candidates = []
    for line in pdb_or_pdbqt.read_text().splitlines():
        if line.startswith("HETATM") or line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            elem = line[76:78].strip() if len(line) > 78 else ""
            # Filter out water O — accept only those with O elem and named O1
            if (elem in ("O", "OA")) and atom_name == "O1":
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                candidates.append(np.array([x, y, z]))
    return candidates[0] if candidates else None


def find_carbonyl_O_in_mode(mode_n: int):
    """Locate the chromenone C=O of α-naphthoflavone in a given mode's SDF
    (SDF preserves bond orders so RDKit can apply ring-aware SMARTS).
    Falls back to the PDB if SDF is unavailable.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    base = VAL / f"redock2_CYP1B1_BHF_mode{mode_n}"
    sdf_path = base.with_suffix(".sdf")
    if not sdf_path.exists():
        # Generate it from the PDBQT mode
        pq = base.with_suffix(".pdbqt")
        if not pq.exists():
            return None
        subprocess.run([OBABEL, str(pq), "-O", str(sdf_path)],
                        capture_output=True, text=True)
    if not sdf_path.exists():
        return None
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next(iter(suppl), None)
    if mol is None:
        return None
    patt = Chem.MolFromSmarts("[O;X1]=[#6;r6]")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        patt = Chem.MolFromSmarts("[O;X1]=[#6]")
        matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None
    o_idx = matches[0][0]
    conf = mol.GetConformer()
    p = conf.GetAtomPosition(o_idx)
    return np.array([p.x, p.y, p.z])


def main() -> int:
    fe = find_heme_iron(RECEPTOR)
    print(f"Heme Fe in 4I8V: {fe}")
    o_xtal = find_chromen4one_carbonyl_O(COCRYS)
    print(f"BHF crystal carbonyl O1: {o_xtal}")
    if fe is None or o_xtal is None:
        print("Could not identify Fe or O1 — aborting")
        return 1
    d_xtal = float(np.linalg.norm(fe - o_xtal))
    print(f"\nCrystal Fe...O(C4=O) distance: {d_xtal:.2f} A\n")

    rows = []
    for mode in range(1, 10):
        mode_pdb = VAL / f"redock2_CYP1B1_BHF_mode{mode}.pdb"
        if not mode_pdb.exists():
            continue
        o_pose = find_carbonyl_O_in_mode(mode)
        if o_pose is None:
            continue
        d = float(np.linalg.norm(fe - o_pose))
        rmsd_r = subprocess.run([OBRMS, str(COCRYS), str(mode_pdb)],
                                capture_output=True, text=True)
        rmsd = None
        try:
            rmsd = float(rmsd_r.stdout.strip().splitlines()[-1].split()[-1])
        except Exception:
            pass
        # Affinity from PDBQT
        aff = None
        pdbqt_mode = VAL / f"redock2_CYP1B1_BHF_mode{mode}.pdbqt"
        if pdbqt_mode.exists():
            for ln in pdbqt_mode.read_text().splitlines():
                if ln.startswith("REMARK VINA RESULT:"):
                    try: aff = float(ln.split()[3]); break
                    except: pass
        within = abs(d - d_xtal) < CONSTRAINT_TOLERANCE_A
        rows.append({
            "mode": mode, "affinity": aff,
            "fe_to_carbonyl_A": round(d, 2),
            "delta_from_xtal": round(d - d_xtal, 2),
            "passes_constraint": bool(within),
            "rmsd_to_crystal": round(rmsd, 2) if rmsd is not None else None,
        })

    df = pd.DataFrame(rows)
    df.to_csv(VAL / "cyp1b1_pharmacophore_filter.tsv", sep="\t", index=False)
    print(df.to_string(index=False))

    constrained = df[df["passes_constraint"]].copy()
    if not constrained.empty:
        # Best-scoring among constraint-passers (most negative affinity)
        best = constrained.loc[constrained["affinity"].idxmin()]
        print(f"\n*** Constraint-satisfying best-scoring mode: #{int(best['mode'])} ***")
        print(f"  affinity:             {best['affinity']}")
        print(f"  Fe...carbonyl_O:      {best['fe_to_carbonyl_A']} A "
              f"(crystal: {d_xtal:.2f} A)")
        print(f"  RMSD to crystal:      {best['rmsd_to_crystal']} A")
    else:
        print("\nNo modes pass the +/- {} A constraint!".format(CONSTRAINT_TOLERANCE_A))

    summary = (
        f"CYP1B1 Pharmacophore-Constrained Redock\n"
        f"========================================\n\n"
        f"Constraint:  Fe ... C4-carbonyl O distance within {CONSTRAINT_TOLERANCE_A} A "
        f"of crystal ({d_xtal:.2f} A)\n\n"
    )
    if not constrained.empty:
        b = constrained.loc[constrained["affinity"].idxmin()]
        summary += f"Mode selected:        #{int(b['mode'])}\n"
        summary += f"Affinity:             {b['affinity']} kcal/mol\n"
        summary += f"Fe...carbonyl_O:      {b['fe_to_carbonyl_A']} A\n"
        summary += f"RMSD vs crystal:      {b['rmsd_to_crystal']} A\n\n"
        summary += f"vs unconstrained mode 1: {df.iloc[0]['rmsd_to_crystal']} A\n"
    (VAL / "cyp1b1_pharmacophore_summary.txt").write_text(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
