"""Run AutoDock Vina across all (ligand, receptor) jobs in dock_jobs.tsv.

For each job:
  vina --receptor R.pdbqt --ligand L.pdbqt
       --center_x ... --size_x ...
       --exhaustiveness 16 --num_modes 9 --seed 42
       --out output/docking/results/{cid}_{pdb}_out.pdbqt
       --log output/docking/results/{cid}_{pdb}.log

Parses each log for the 9 modes' affinity (kcal/mol) and rmsd_lb / rmsd_ub,
then summarises in dock_results.tsv with the top-1 affinity per (compound, target).
"""

from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
RESULTS = DOCK / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

VINA = "/usr/local/bin/vina"
EXHAUSTIVENESS = 16
NUM_MODES = 9


def run_one(job: dict) -> dict:
    cid = job["compound_id"]
    pdb = job["pdb"]
    out_pdbqt = RESULTS / f"{cid}_{pdb}_out.pdbqt"
    log = RESULTS / f"{cid}_{pdb}.log"
    if log.exists() and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
        return parse_log(log, job)
    cmd = [
        VINA,
        "--receptor", job["receptor_pdbqt"],
        "--ligand",   job["ligand_pdbqt"],
        "--center_x", str(job["center_x"]),
        "--center_y", str(job["center_y"]),
        "--center_z", str(job["center_z"]),
        "--size_x",   str(job["size_x"]),
        "--size_y",   str(job["size_y"]),
        "--size_z",   str(job["size_z"]),
        "--out",      str(out_pdbqt),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", "42",
    ]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    log.write_text(r.stdout + "\n" + r.stderr)
    if r.returncode != 0:
        print(f"  FAIL {cid} {pdb}: rc={r.returncode}", file=sys.stderr)
        print("  STDERR:", r.stderr[:400], file=sys.stderr)
        return {**job, "best_affinity": None, "elapsed_s": round(elapsed, 1), "n_modes": 0}
    parsed = parse_log(log, job)
    parsed["elapsed_s"] = round(elapsed, 1)
    return parsed


MODE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$")


def parse_log(log_path: Path, job: dict) -> dict:
    modes = []
    in_table = False
    for line in log_path.read_text().splitlines():
        if "mode |   affinity" in line.lower() or "mode |  affinity" in line.lower():
            in_table = True
            continue
        if in_table:
            m = MODE_RE.match(line)
            if m:
                idx, aff, rl, ru = m.groups()
                modes.append({
                    "mode": int(idx),
                    "affinity": float(aff),
                    "rmsd_lb": float(rl),
                    "rmsd_ub": float(ru),
                })
            elif line.strip().startswith("Writing") or line.strip() == "":
                if modes:
                    in_table = False
    best = min((m["affinity"] for m in modes), default=None)
    return {
        "compound_id": job["compound_id"],
        "target": job["target"],
        "uniprot": job["uniprot"],
        "pdb": job["pdb"],
        "n_modes": len(modes),
        "best_affinity": best,
        "mode2_affinity": modes[1]["affinity"] if len(modes) > 1 else None,
        "mode3_affinity": modes[2]["affinity"] if len(modes) > 2 else None,
    }


def main() -> int:
    jobs = pd.read_csv(DOCK / "dock_jobs.tsv", sep="\t").to_dict("records")
    print(f"Running {len(jobs)} dock jobs ...", file=sys.stderr)
    rows = []
    t_start = time.time()
    for i, job in enumerate(jobs, 1):
        out = run_one(job)
        rows.append(out)
        elapsed = time.time() - t_start
        eta = elapsed / i * (len(jobs) - i)
        print(f"  [{i:>3}/{len(jobs)}] {out['compound_id']:<10} {out['target']:<8} "
              f"best={out.get('best_affinity', 'NA')}  modes={out['n_modes']}  "
              f"elapsed={elapsed:.0f}s eta={eta:.0f}s",
              file=sys.stderr)
    df = pd.DataFrame(rows).sort_values(["target", "best_affinity"])
    df.to_csv(DOCK / "dock_results.tsv", sep="\t", index=False)
    print(f"\nWrote {DOCK / 'dock_results.tsv'} ({len(df)} rows)", file=sys.stderr)

    # Summary
    print("\n=== Summary by target ===", file=sys.stderr)
    for tname, grp in df.groupby("target"):
        best = grp.loc[grp["best_affinity"].idxmin()]
        print(f"  {tname}  n={len(grp)}  median={grp['best_affinity'].median():.2f} "
              f"top={best['best_affinity']:.2f} ({best['compound_id']})",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
