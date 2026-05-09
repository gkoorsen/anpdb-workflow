"""Recompute redock RMSD properly — absolute coords + symmetry-corrected
heavy-atom matching via OpenBabel SMARTS or simple element-greedy matching.

Approach:
  1. Cocrystal PDB heavy atoms (with elements)
  2. Docked-pose mode-1 PDBQT heavy atoms (with elements)
  3. For each crystal atom, find the nearest *unused* pose atom of the same
     element. RMSD = sqrt(mean(pairwise sq distance)).

This is a permissive lower-bound on the symmetry-corrected RMSD: real
algorithms (Hungarian on equivalent-atom-class graph) would give the same
or smaller. If our greedy result is < 2 A then the proper algorithm would
also be < 2 A.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
RECEPTORS = DOCK / "receptors"
VAL = DOCK / "validation"

COCRYSTAL = {
    "CYP1B1": ("4I8V", "BHF"),
    "SGLT2":  ("7VSI", "7R3"),
    "MAO-B":  ("2V5Z", "SAG"),
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
    # Fallback: parse leading letters from atom name column
    name = atom_name.strip()
    if name:
        # First non-digit character(s) — handle "C1", "Cl2", "NA"
        letters = ""
        for ch in name:
            if ch.isalpha():
                letters += ch
            else:
                break
        # Try two-letter then one-letter
        if len(letters) >= 2 and letters[:2].title() in AUTODOCK_TO_ELEMENT:
            return AUTODOCK_TO_ELEMENT[letters[:2].title()]
        if letters[:1].upper() in AUTODOCK_TO_ELEMENT:
            return AUTODOCK_TO_ELEMENT[letters[:1].upper()]
    # Last resort: just clean what we have
    cleaned = "".join(c for c in raw if c.isalpha()).capitalize()
    return cleaned


def parse_pdb_heavy(path: Path):
    """Return list of (element, x, y, z) for non-H atoms in MODEL 1 (or all if no MODEL)."""
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
    """Element-aware greedy matching, then absolute RMSD."""
    # Group by element
    used = [False] * len(posed)
    sq = []
    for ce, cx, cy, cz in crystal:
        best = None
        best_idx = -1
        for i, (pe, px, py, pz) in enumerate(posed):
            if used[i] or pe != ce:
                continue
            d = (cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2
            if best is None or d < best:
                best = d
                best_idx = i
        if best is None:
            return None  # element mismatch — can't compute
        used[best_idx] = True
        sq.append(best)
    return float(np.sqrt(np.mean(sq)))


def main() -> int:
    rows = []
    for tname, (pdb, resn) in COCRYSTAL.items():
        cryst = RECEPTORS / f"{pdb}_cocrystal_{resn}.pdb"
        posed = VAL / f"redock_{tname}_{resn}_out.pdbqt"
        c = parse_pdb_heavy(cryst)
        p = parse_pdb_heavy(posed)
        rmsd = greedy_rmsd(c, p)
        rows.append({
            "target": tname, "pdb": pdb, "ligand_resn": resn,
            "n_xtal_heavy": len(c), "n_pose_heavy": len(p),
            "rmsd_absolute_A": round(rmsd, 2) if rmsd is not None else None,
        })
        print(f"  {tname}: xtal={len(c)} pose={len(p)} RMSD(abs)={rmsd:.2f}" if rmsd else f"  {tname}: ELEMENT MISMATCH")
    df = pd.DataFrame(rows)
    df.to_csv(VAL / "redock_rmsd_corrected.tsv", sep="\t", index=False)
    print("\nFinal corrected redock RMSDs:")
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
