"""Combined MD analysis figure: backbone RMSD + ligand-binding-site RMSF
for both CYP1B1 and SGLT2 production trajectories.

Usage: python md_analysis.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
MD_ROOT = ROOT / "output" / "md"

sns.set_theme(context="talk", style="whitegrid", palette="deep")

TARGET_COLOURS = {"CYP1B1": "#1C7293", "SGLT2": "#EE854A"}


def analyse_target(target: str, top_pdb: Path, traj_dcd: Path, lead_id: str):
    """Compute backbone RMSD + ligand-pocket residue RMSF.

    Returns a dict ready for plotting.
    """
    import mdtraj as md

    traj = md.load(str(traj_dcd), top=str(top_pdb))
    print(f"[{target}] frames={traj.n_frames}  atoms={traj.n_atoms}",
          file=sys.stderr)

    # Strip waters/ions for analysis
    prot_lig = traj.atom_slice(traj.topology.select("not (resname HOH or resname NA or resname CL)"))

    # Backbone RMSD (against frame 0, no realignment of ligand)
    backbone = prot_lig.topology.select("backbone")
    prot_lig.superpose(prot_lig, frame=0, atom_indices=backbone)
    rmsd = md.rmsd(prot_lig, prot_lig, frame=0, atom_indices=backbone) * 10  # nm -> Å

    # Ligand atoms — anything not protein/water/ion
    lig_idx = prot_lig.topology.select(
        "not protein and not (resname HOH or resname NA or resname CL)"
    )
    if len(lig_idx) == 0:
        # Fallback: explicit residue name UNL or LIG
        lig_idx = prot_lig.topology.select("resname UNL or resname LIG or resname UNK")
    print(f"[{target}] ligand atoms in trajectory: {len(lig_idx)}", file=sys.stderr)

    # Ligand RMSD vs frame 0 (after backbone alignment to frame 0)
    if len(lig_idx) > 0:
        lig_rmsd = md.rmsd(prot_lig, prot_lig, frame=0, atom_indices=lig_idx) * 10
    else:
        lig_rmsd = None

    # Residue RMSF (Calpha)
    ca = prot_lig.topology.select("name CA")
    rmsf = md.rmsf(prot_lig, prot_lig, frame=0, atom_indices=ca) * 10  # nm -> Å

    times = np.arange(traj.n_frames) * 0.01  # 5000-step interval × 0.002 ps/step = 10 ps/frame = 0.01 ns

    return {
        "target": target,
        "lead": lead_id,
        "n_frames": traj.n_frames,
        "times_ns": times,
        "backbone_rmsd": rmsd,
        "ligand_rmsd": lig_rmsd,
        "rmsf_ca": rmsf,
        "ca_residue_indices": np.arange(len(ca)),
    }


def main():
    targets = []
    for tname in ["CYP1B1", "SGLT2"]:
        d = MD_ROOT / tname
        traj = d / "md_trajectory.dcd"
        top = d / "md_solvated.pdb"
        # Find the lead compound id from a pose file
        lead_pdbs = list(d.glob("Mol_*_pose.pdb"))
        if not traj.exists() or not top.exists() or not lead_pdbs:
            print(f"[{tname}] missing trajectory/topology — skipping", file=sys.stderr)
            continue
        lead = lead_pdbs[0].stem.replace("_pose", "")
        result = analyse_target(tname, top, traj, lead)
        targets.append(result)

    if not targets:
        print("No completed MD trajectories", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))

    # Top-left: backbone RMSD over time, both targets
    ax = axes[0, 0]
    for t in targets:
        ax.plot(t["times_ns"], t["backbone_rmsd"],
                label=f"{t['target']} ({t['lead']})",
                color=TARGET_COLOURS[t["target"]], lw=1.5)
    ax.axhline(2.5, color="gray", ls="--", lw=1, alpha=0.5)
    ax.text(t["times_ns"][-1] * 0.98, 2.55, "2.5 Å (typical stable)",
            ha="right", fontsize=9, color="gray", italic=True if False else False)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Backbone RMSD (Å)")
    ax.set_title("(A) Protein backbone stability", loc="left", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=10)

    # Top-right: ligand RMSD over time
    ax = axes[0, 1]
    for t in targets:
        if t["ligand_rmsd"] is not None:
            ax.plot(t["times_ns"], t["ligand_rmsd"],
                    label=f"{t['target']} ({t['lead']})",
                    color=TARGET_COLOURS[t["target"]], lw=1.5)
    ax.axhline(3.0, color="gray", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("(B) Ligand pose stability (vs frame 0)", loc="left",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=10)

    # Bottom-left: CYP1B1 RMSF if present
    cyp = next((t for t in targets if t["target"] == "CYP1B1"), None)
    sglt = next((t for t in targets if t["target"] == "SGLT2"), None)

    ax = axes[1, 0]
    if cyp is not None:
        ax.plot(cyp["ca_residue_indices"], cyp["rmsf_ca"],
                color=TARGET_COLOURS["CYP1B1"], lw=1)
        ax.fill_between(cyp["ca_residue_indices"], 0, cyp["rmsf_ca"],
                        color=TARGET_COLOURS["CYP1B1"], alpha=0.3)
        ax.set_title(f"(C) CYP1B1 Cα RMSF ({cyp['lead']})", loc="left",
                     fontsize=12, fontweight="bold")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("RMSF (Å)")

    ax = axes[1, 1]
    if sglt is not None:
        ax.plot(sglt["ca_residue_indices"], sglt["rmsf_ca"],
                color=TARGET_COLOURS["SGLT2"], lw=1)
        ax.fill_between(sglt["ca_residue_indices"], 0, sglt["rmsf_ca"],
                        color=TARGET_COLOURS["SGLT2"], alpha=0.3)
        ax.set_title(f"(D) SGLT2 Cα RMSF ({sglt['lead']})", loc="left",
                     fontsize=12, fontweight="bold")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("RMSF (Å)")

    fig.suptitle(
        "Molecular dynamics validation — top consensus leads\n"
        "2 ns NPT @ 310 K, AMBER ff14SB + SMIRNOFF (NAGL), TIP3P, 0.15 M NaCl",
        fontsize=12, fontweight="bold", y=0.99,
    )
    fig.tight_layout()

    out_path = MD_ROOT / "fig_md_combined.png"
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)

    # Save numeric data
    summary = "MD Analysis Summary\n"
    summary += "=" * 40 + "\n\n"
    for t in targets:
        summary += f"{t['target']} ({t['lead']})\n"
        summary += f"  Frames:                {t['n_frames']}\n"
        summary += f"  Backbone RMSD final:   {t['backbone_rmsd'][-1]:.2f} Å\n"
        summary += f"  Backbone RMSD mean:    {t['backbone_rmsd'].mean():.2f} Å\n"
        summary += f"  Backbone RMSD max:     {t['backbone_rmsd'].max():.2f} Å\n"
        if t["ligand_rmsd"] is not None:
            summary += f"  Ligand RMSD final:     {t['ligand_rmsd'][-1]:.2f} Å\n"
            summary += f"  Ligand RMSD mean:      {t['ligand_rmsd'].mean():.2f} Å\n"
            summary += f"  Ligand RMSD max:       {t['ligand_rmsd'].max():.2f} Å\n"
        summary += f"  Cα RMSF mean:          {t['rmsf_ca'].mean():.2f} Å\n"
        summary += f"  Cα RMSF max:           {t['rmsf_ca'].max():.2f} Å\n\n"

    (MD_ROOT / "md_combined_summary.txt").write_text(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
