"""Best-of-top-9-modes RMSD analysis for all 3 redocks (SGLT2 + MAO-B + CYP1B1).

Reports per-mode affinity + RMSD vs crystal, and the "best-of-9" RMSD.
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

REDOCKS = [
    ("CYP1B1", "4I8V", "BHF"),
    ("SGLT2",  "7VSI", "7R3"),
    ("MAO-B",  "2V5Z", "SAG"),
]


def split_modes(pdbqt: Path, prefix: str):
    text = pdbqt.read_text().splitlines()
    paths = []
    cur_lines = None
    cur_n = None
    cur_aff = None
    for line in text:
        if line.startswith("MODEL"):
            cur_n = int(line.split()[1])
            cur_aff = None
            cur_lines = [line]
        elif line.startswith("REMARK VINA RESULT:") and cur_lines is not None:
            try: cur_aff = float(line.split()[3])
            except: pass
            cur_lines.append(line)
        elif line.startswith("ENDMDL") and cur_lines is not None:
            cur_lines.append(line)
            tmp_pq = VAL / f"{prefix}_mode{cur_n}.pdbqt"
            tmp_pdb = VAL / f"{prefix}_mode{cur_n}.pdb"
            tmp_pq.write_text("\n".join(cur_lines))
            subprocess.run([OBABEL, str(tmp_pq), "-O", str(tmp_pdb)],
                           capture_output=True)
            paths.append((cur_n, cur_aff, tmp_pdb))
            cur_lines = None
        elif cur_lines is not None:
            cur_lines.append(line)
    return paths


def rmsd_obrms(reference: Path, pose: Path):
    r = subprocess.run([OBRMS, str(reference), str(pose)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    line = r.stdout.strip().splitlines()[-1]
    try:
        return float(line.split()[-1])
    except (ValueError, IndexError):
        return None


def main():
    rows_per_mode = []
    summary_rows = []
    for tname, pdb, resn in REDOCKS:
        out = VAL / f"redock2_{tname}_{resn}_out.pdbqt"
        if not out.exists():
            print(f"  SKIP {tname}: no PDBQT")
            continue
        cocrys = RECEPTORS / f"{pdb}_cocrystal_{resn}.pdb"
        modes = split_modes(out, prefix=f"redock2_{tname}_{resn}")
        print(f"\n=== {tname} ({pdb}/{resn}) — {len(modes)} modes ===")
        per_target = []
        for n, aff, mode_pdb in modes:
            rmsd = rmsd_obrms(cocrys, mode_pdb)
            rows_per_mode.append({
                "target": tname, "pdb": pdb, "ligand_resn": resn,
                "mode": n, "affinity": aff,
                "rmsd_to_crystal": round(rmsd, 2) if rmsd is not None else None,
            })
            per_target.append((n, aff, rmsd))
            print(f"  mode {n}: aff={aff:>6}  RMSD={rmsd:>6.2f}" if rmsd is not None
                  else f"  mode {n}: aff={aff:>6}  RMSD=NA")
        valid = [(n, a, r) for n, a, r in per_target if r is not None]
        if not valid:
            continue
        best_score = valid[0]   # mode 1 is always best-scoring
        best_rmsd = min(valid, key=lambda x: x[2])
        summary_rows.append({
            "target": tname,
            "best_score_mode": best_score[0],
            "best_score_affinity": best_score[1],
            "best_score_rmsd": round(best_score[2], 2),
            "best_rmsd_mode": best_rmsd[0],
            "best_rmsd_affinity": best_rmsd[1],
            "best_rmsd": round(best_rmsd[2], 2),
            "score_gap_kcal": round(best_score[1] - best_rmsd[1], 2),
            "n_modes_under_2A": sum(1 for _, _, r in valid if r < 2.0),
            "n_modes_under_3A": sum(1 for _, _, r in valid if r < 3.0),
        })
        print(f"  Best-scoring mode #{best_score[0]}: "
              f"aff={best_score[1]} RMSD={best_score[2]:.2f}")
        print(f"  Best-RMSD mode    #{best_rmsd[0]}: "
              f"aff={best_rmsd[1]} RMSD={best_rmsd[2]:.2f}")
        print(f"  Modes <2A: {summary_rows[-1]['n_modes_under_2A']}/{len(valid)}")

    pd.DataFrame(rows_per_mode).to_csv(VAL / "redock_all_modes.tsv", sep="\t", index=False)
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(VAL / "redock_best_of_9.tsv", sep="\t", index=False)
    print("\n\nSummary table:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
