"""Generate publication figures for the docking validation results.

Reads output/docking/dock_results.tsv and produces:
- fig_dock_boxplot.png         per-target affinity boxplot
- fig_dock_top_leads.png       per-target leaderboard (top hits as bars)
- fig_dock_summary.png         combined: stripplot + medians + receptor inset
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
DOCK = ROOT / "output" / "docking"

sns.set_theme(context="talk", style="whitegrid", palette="deep")

TARGET_COLOURS = {"CYP1B1": "#1C7293", "SGLT2": "#EE854A", "MAO-B": "#2C8C5A"}


def main() -> int:
    df = pd.read_csv(DOCK / "dock_results.tsv", sep="\t")
    df = df.dropna(subset=["best_affinity"]).copy()
    # Drop the two failed runs (n_modes <= 1) — keep meaningful poses only
    df = df[df["n_modes"] >= 4].copy()
    print(f"Plotting {len(df)} valid docking results", end=" ")
    print(f"({df['target'].value_counts().to_dict()})")

    target_order = ["CYP1B1", "SGLT2", "MAO-B"]

    # --- Strip + box per target ---
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.boxplot(
        data=df, x="target", y="best_affinity",
        order=target_order, palette=TARGET_COLOURS,
        width=0.5, fliersize=0, ax=ax,
        boxprops={"alpha": 0.4},
    )
    sns.stripplot(
        data=df, x="target", y="best_affinity",
        order=target_order, palette=TARGET_COLOURS,
        size=8, edgecolor="white", linewidth=0.7, ax=ax,
        jitter=0.18,
    )
    # Annotate top hits and counts
    for i, t in enumerate(target_order):
        sub = df[df["target"] == t]
        if sub.empty:
            continue
        med = sub["best_affinity"].median()
        n = len(sub)
        best = sub["best_affinity"].min()
        best_cid = sub.loc[sub["best_affinity"].idxmin(), "compound_id"]
        ax.text(i, ax.get_ylim()[1] - 0.5,
                f"n={n}\nmedian={med:.1f}",
                ha="center", fontsize=10, color=TARGET_COLOURS[t], fontweight="bold")
        ax.scatter([i], [best], marker="*", s=300, color="#D62728",
                   edgecolor="white", linewidth=1.5, zorder=10)
        ax.annotate(f"{best_cid}\n({best:.1f})", (i, best),
                    xytext=(15, -10), textcoords="offset points",
                    fontsize=9, color="#D62728", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#D62728", lw=1))

    ax.axhline(-7, color="gray", lw=1, ls="--", alpha=0.5)
    ax.text(2.4, -7.1, "−7 kcal/mol\n(typical hit threshold)",
            ha="right", fontsize=9, color="gray", style="italic")
    ax.set_ylabel("Best Vina affinity (kcal/mol)\n← stronger binding")
    ax.set_xlabel("")
    ax.set_title(
        "Docking validation of 33 consensus leads against 3 targets\n"
        "(AutoDock Vina, exhaustiveness=16, 9 modes per ligand)",
        loc="left", fontsize=13,
    )
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(DOCK / "fig_dock_boxplot.png", dpi=160)
    plt.close(fig)
    print(f"  -> fig_dock_boxplot.png")

    # --- Per-target top-N bars ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=False)
    for ax, t in zip(axes, target_order):
        sub = df[df["target"] == t].sort_values("best_affinity").head(8)
        if sub.empty:
            ax.set_visible(False)
            continue
        y = np.arange(len(sub))
        ax.barh(y, sub["best_affinity"], color=TARGET_COLOURS[t], alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["compound_id"], fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("Best Vina affinity (kcal/mol)")
        ax.set_title(f"{t} ({sub.iloc[0]['pdb']})\nTop {len(sub)} of {(df['target']==t).sum()}",
                     loc="left", fontsize=12, color=TARGET_COLOURS[t], fontweight="bold")
        for yi, val in enumerate(sub["best_affinity"]):
            ax.text(val - 0.15, yi, f"{val:.1f}",
                    va="center", ha="right", fontsize=9, color="white", fontweight="bold")
        ax.set_xlim(min(sub["best_affinity"].min() - 0.5, -8), 0)
    fig.suptitle("Per-target docking leaderboards", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(DOCK / "fig_dock_top_leads.png", dpi=160)
    plt.close(fig)
    print(f"  -> fig_dock_top_leads.png")

    # --- Summary: print the 3 best hits ---
    print("\nTop hit per target:")
    for t in target_order:
        sub = df[df["target"] == t]
        if sub.empty:
            continue
        best = sub.loc[sub["best_affinity"].idxmin()]
        print(f"  {t} ({best['pdb']}): {best['compound_id']}  "
              f"affinity={best['best_affinity']:.2f} kcal/mol  "
              f"(median across {len(sub)} cpds = {sub['best_affinity'].median():.2f})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
