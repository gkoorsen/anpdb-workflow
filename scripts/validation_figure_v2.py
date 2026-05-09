"""Updated validation figure including best-of-9 RMSD analysis.

Three-panel figure:
  (A) Decoy control: consensus vs property-matched decoys per target
  (B) Best-of-top-9 RMSD recovery for all 3 redocks
  (C) CYP1B1 mode-by-mode scatter — affinity vs RMSD, highlighting the
      best-RMSD pose (mode 6) and the best-scoring pose (mode 1)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"
VAL = DOCK / "validation"

sns.set_theme(context="talk", style="whitegrid", palette="deep")

TARGET_COLOURS = {"CYP1B1": "#1C7293", "SGLT2": "#EE854A", "MAO-B": "#2C8C5A"}


def main():
    consensus = pd.read_csv(DOCK / "dock_results.tsv", sep="\t").dropna(subset=["best_affinity"])
    consensus = consensus[consensus["n_modes"] >= 4]
    consensus["set"] = "Consensus hit"

    decoy = pd.read_csv(VAL / "decoy_results.tsv", sep="\t").dropna(subset=["best_affinity"])
    decoy["set"] = "Property-matched decoy"

    combined = pd.concat([
        consensus[["target", "best_affinity", "set"]],
        decoy[["target", "best_affinity", "set"]],
    ], ignore_index=True)

    redock = pd.read_csv(VAL / "redock_best_of_9.tsv", sep="\t")
    cyp_modes = pd.read_csv(VAL / "redock_all_modes.tsv", sep="\t")
    cyp_modes = cyp_modes[cyp_modes["target"] == "CYP1B1"]

    target_order = ["CYP1B1", "SGLT2", "MAO-B"]

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 6, hspace=0.45, wspace=0.5)

    # (A) Decoy control
    ax0 = fig.add_subplot(gs[0, 0:3])
    sns.boxplot(
        data=combined, x="target", y="best_affinity", hue="set",
        order=target_order, hue_order=["Consensus hit", "Property-matched decoy"],
        palette={"Consensus hit": "#1C7293", "Property-matched decoy": "#AAAAAA"},
        width=0.55, ax=ax0, fliersize=0,
        boxprops={"alpha": 0.5},
    )
    sns.stripplot(
        data=combined, x="target", y="best_affinity", hue="set",
        order=target_order, hue_order=["Consensus hit", "Property-matched decoy"],
        palette={"Consensus hit": "#1C7293", "Property-matched decoy": "#666666"},
        size=5, edgecolor="white", linewidth=0.4, dodge=True, ax=ax0, jitter=0.18,
        legend=False,
    )
    ax0.axhline(-7, color="gray", lw=1, ls="--", alpha=0.5)
    ax0.set_ylabel("Best Vina affinity\n(kcal/mol)")
    ax0.set_xlabel("")
    ax0.set_title("(A) Decoy control", loc="left", fontsize=13, fontweight="bold")
    ax0.invert_yaxis()
    ax0.legend(loc="upper right", title="", fontsize=10)

    # (B) Best-of-9 RMSD comparison
    ax1 = fig.add_subplot(gs[0, 3:])
    x = np.arange(len(redock))
    w = 0.4
    bs = ax1.bar(x - w/2, redock["best_score_rmsd"], w,
                 label="Best-scoring (mode 1)",
                 color=[TARGET_COLOURS[t] for t in redock["target"]], alpha=0.5)
    br = ax1.bar(x + w/2, redock["best_rmsd"], w,
                 label="Best-of-top-9",
                 color=[TARGET_COLOURS[t] for t in redock["target"]], alpha=1.0)
    ax1.axhline(2.0, color="#D62728", lw=1.5, ls="--", alpha=0.7,
                label="2 Å threshold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(redock["target"])
    ax1.set_ylabel("RMSD vs crystal pose (Å)")
    for i, (_, r) in enumerate(redock.iterrows()):
        ax1.text(i - w/2, r["best_score_rmsd"] + 0.1,
                 f"{r['best_score_rmsd']:.2f}",
                 ha="center", fontsize=9, color="#444")
        ax1.text(i + w/2, r["best_rmsd"] + 0.1,
                 f"{r['best_rmsd']:.2f}",
                 ha="center", fontsize=10, fontweight="bold", color="#222")
    ax1.set_title("(B) Redock control — best-of-top-9 within 2 Å for all 3 targets",
                  loc="left", fontsize=13, fontweight="bold")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.set_ylim(0, max(redock["best_score_rmsd"]) + 1.5)

    # (C) CYP1B1 mode-by-mode scatter
    ax2 = fig.add_subplot(gs[1, :])
    ax2.scatter(cyp_modes["rmsd_to_crystal"], -cyp_modes["affinity"],
                s=200, c=TARGET_COLOURS["CYP1B1"], edgecolor="white", linewidth=1.2, alpha=0.8)
    for _, r in cyp_modes.iterrows():
        ax2.annotate(f"#{int(r['mode'])}",
                     (r["rmsd_to_crystal"], -r["affinity"]),
                     xytext=(4, 4), textcoords="offset points",
                     fontsize=11, fontweight="bold")

    # Highlight mode 1 (best-scoring) and best-RMSD mode
    best_score = cyp_modes.iloc[cyp_modes["affinity"].idxmin()]
    best_rmsd  = cyp_modes.iloc[cyp_modes["rmsd_to_crystal"].idxmin()]
    ax2.scatter([best_score["rmsd_to_crystal"]], [-best_score["affinity"]],
                s=400, marker="*", color="#D62728", edgecolor="white", linewidth=1.5,
                zorder=10, label=f"Best-scoring: mode {int(best_score['mode'])} ({best_score['rmsd_to_crystal']:.2f} Å, flipped)")
    ax2.scatter([best_rmsd["rmsd_to_crystal"]], [-best_rmsd["affinity"]],
                s=400, marker="*", color="#2CA02C", edgecolor="white", linewidth=1.5,
                zorder=10, label=f"Best-RMSD: mode {int(best_rmsd['mode'])} ({best_rmsd['rmsd_to_crystal']:.2f} Å, crystal pose)")

    ax2.axvline(2.0, color="#D62728", lw=1, ls="--", alpha=0.5)
    ax2.text(2.1, ax2.get_ylim()[1] - 0.5, "2 Å",
             color="#D62728", fontsize=10)
    ax2.set_xlabel("RMSD vs crystal pose (Å, symmetry-corrected via obrms)")
    ax2.set_ylabel("|Vina affinity| (kcal/mol)\n← stronger binding")
    ax2.set_title(
        "(C) CYP1B1 / α-naphthoflavone — all 9 Vina modes\n"
        "Crystal pose (mode 6, RMSD 0.39 Å) ranks 4th by score; flipped pose (mode 1) "
        "scores 1.0 kcal/mol better — symmetric-flip artefact for flat aromatic ligands.",
        loc="left", fontsize=12, fontweight="bold",
    )
    ax2.legend(loc="lower right", fontsize=10)

    fig.suptitle(
        "Docking validation controls",
        fontsize=15, fontweight="bold", y=0.99,
    )
    fig.tight_layout()
    fig.savefig(VAL / "fig_dock_validation_v2.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {VAL / 'fig_dock_validation_v2.png'}")


if __name__ == "__main__":
    main()
