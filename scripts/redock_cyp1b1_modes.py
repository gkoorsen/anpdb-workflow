"""Check ALL 9 Vina modes for the CYP1B1 redock — was the crystal pose
sampled but mis-ranked, or genuinely missed?

The standard literature solution to a symmetric ring-flip artifact in
docking validation is to report the BEST-RMSD pose across the top-N
modes, not the best-scoring mode. If the protocol samples the crystal
pose (RMSD < 2 A) anywhere in modes 1–9, the issue is a scoring-function
ranking problem, not a sampling problem.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
VAL = ROOT / "output" / "docking" / "validation"
RECEPTORS = ROOT / "output" / "docking" / "receptors"
OBABEL = "/opt/homebrew/bin/obabel"
OBRMS = "/opt/homebrew/bin/obrms"

REDOCK_PDBQT = VAL / "redock2_CYP1B1_BHF_out.pdbqt"
COCRYSTAL = RECEPTORS / "4I8V_cocrystal_BHF.pdb"


def split_modes_to_pdb():
    """Split each MODEL of the multi-pose PDBQT into a separate PDB."""
    text = REDOCK_PDBQT.read_text().splitlines()
    paths = []
    current_lines = None
    current_n = None
    for line in text:
        if line.startswith("MODEL"):
            current_n = int(line.split()[1])
            current_lines = [line]
        elif line.startswith("ENDMDL"):
            current_lines.append(line)
            tmp_pdbqt = VAL / f"redock2_CYP1B1_mode{current_n}.pdbqt"
            tmp_pdb = VAL / f"redock2_CYP1B1_mode{current_n}.pdb"
            tmp_pdbqt.write_text("\n".join(current_lines))
            subprocess.run([OBABEL, str(tmp_pdbqt), "-O", str(tmp_pdb)],
                           capture_output=True)
            paths.append((current_n, tmp_pdb))
            current_lines = None
        elif current_lines is not None:
            current_lines.append(line)
    return paths


def affinity_per_mode():
    """Parse REMARK VINA RESULT per MODEL."""
    affs = {}
    text = REDOCK_PDBQT.read_text().splitlines()
    current_n = None
    for line in text:
        if line.startswith("MODEL"):
            current_n = int(line.split()[1])
        elif line.startswith("REMARK VINA RESULT:") and current_n is not None:
            try:
                affs[current_n] = float(line.split()[3])
            except ValueError:
                pass
    return affs


def compute_rmsd(mode_pdb: Path) -> float:
    r = subprocess.run([OBRMS, str(COCRYSTAL), str(mode_pdb)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    # Parse last token of last line
    line = r.stdout.strip().splitlines()[-1]
    try:
        return float(line.split()[-1])
    except (ValueError, IndexError):
        return None


def main() -> int:
    print("Splitting CYP1B1 redock PDBQT into 9 modes ...")
    paths = split_modes_to_pdb()
    print(f"  {len(paths)} modes extracted")

    affs = affinity_per_mode()
    rows = []
    for n, pdb in paths:
        rmsd = compute_rmsd(pdb)
        rows.append({
            "mode": n,
            "affinity": affs.get(n),
            "rmsd_to_crystal": round(rmsd, 2) if rmsd is not None else None,
        })

    df = pd.DataFrame(rows).sort_values("mode")
    df.to_csv(VAL / "redock_CYP1B1_all_modes.tsv", sep="\t", index=False)

    print("\nMode-by-mode RMSD vs crystal pose:")
    print(df.to_string(index=False))

    valid = df.dropna(subset=["rmsd_to_crystal"])
    if not valid.empty:
        best_mode = valid.loc[valid["rmsd_to_crystal"].idxmin()]
        print(f"\nBest-RMSD mode: #{int(best_mode['mode'])}  "
              f"affinity={best_mode['affinity']}  RMSD={best_mode['rmsd_to_crystal']} A")
        n_under_2 = (valid["rmsd_to_crystal"] < 2.0).sum()
        n_under_3 = (valid["rmsd_to_crystal"] < 3.0).sum()
        print(f"\nModes with RMSD < 2 A: {n_under_2}/{len(valid)}")
        print(f"Modes with RMSD < 3 A: {n_under_3}/{len(valid)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
