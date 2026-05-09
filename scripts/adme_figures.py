"""Generate ADME publication figures for the top-50 consensus compounds.

Reads output/adme/adme_top50.csv and produces:
- fig_boiled_egg.png     TPSA vs LogP scatter (BOILED-Egg regions)
- fig_adme_rules.png     pass/fail bar chart for Lipinski/Veber/Egan/GI/BBB
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "adme"

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def main() -> int:
    df = pd.read_csv(OUT / "adme_top50.csv")
    n = len(df)

    # --- BOILED-Egg plot ---
    fig, ax = plt.subplots(figsize=(9, 7))
    # GI absorption region
    gi_rect = mpatches.FancyBboxPatch(
        (-0.4, 0), 6.0, 142, boxstyle="round,pad=0.2",
        facecolor="#FFE0B2", edgecolor="#FF9800", alpha=0.3, lw=1.5,
    )
    ax.add_patch(gi_rect)
    # BBB region
    bbb_rect = mpatches.FancyBboxPatch(
        (-1.0, 0), 4.5, 79, boxstyle="round,pad=0.2",
        facecolor="#BBDEFB", edgecolor="#1976D2", alpha=0.3, lw=1.5,
    )
    ax.add_patch(bbb_rect)

    colors = []
    for _, r in df.iterrows():
        if r["BBB_permeant"]:
            colors.append("#1976D2")
        elif r["GI_absorption"]:
            colors.append("#FF9800")
        else:
            colors.append("#888888")

    ax.scatter(df["LogP"], df["TPSA"], c=colors, s=50, alpha=0.7, edgecolors="white", lw=0.5)

    ax.set_xlabel("WLogP", fontsize=12)
    ax.set_ylabel("TPSA (A²)", fontsize=12)
    ax.set_title("BOILED-Egg: GI absorption and BBB permeability\n(top 50 consensus compounds)",
                 loc="left", fontsize=13)

    legend_elements = [
        mpatches.Patch(facecolor="#BBDEFB", edgecolor="#1976D2", label=f"BBB permeant ({df['BBB_permeant'].sum()}/{n})"),
        mpatches.Patch(facecolor="#FFE0B2", edgecolor="#FF9800", label=f"GI absorbed ({df['GI_absorption'].sum()}/{n})"),
        mpatches.Patch(facecolor="#DDDDDD", edgecolor="#888888", label=f"Low absorption ({n - df['GI_absorption'].sum()}/{n})"),
    ]
    ax.legend(handles=legend_elements, fontsize=10, loc="upper right")
    ax.set_xlim(-7, 10)
    ax.set_ylim(-10, 400)
    fig.tight_layout()
    fig.savefig(OUT / "fig_boiled_egg.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_boiled_egg.png'}")

    # --- Rule pass/fail bar chart ---
    rules = ["Lipinski_pass", "Veber_pass", "Egan_pass", "GI_absorption", "BBB_permeant"]
    labels = ["Lipinski Ro5\n(<=1 violation)", "Veber\n(RotB/TPSA)", "Egan\n(TPSA/LogP)",
              "GI absorption\n(BOILED-Egg)", "BBB permeant\n(BOILED-Egg)"]
    pass_counts = [df[r].sum() for r in rules]
    fail_counts = [n - p for p in pass_counts]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(rules))
    w = 0.35
    ax.bar(x - w/2, pass_counts, w, label="Pass", color="#2CA02C")
    ax.bar(x + w/2, fail_counts, w, label="Fail", color="#D62728", alpha=0.7)
    for i, (p, f) in enumerate(zip(pass_counts, fail_counts)):
        ax.text(i - w/2, p + 0.5, str(p), ha="center", fontsize=10, fontweight="bold", color="#2CA02C")
        ax.text(i + w/2, f + 0.5, str(f), ha="center", fontsize=10, color="#D62728")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Number of compounds")
    ax.set_title(f"ADME rule compliance — top 50 consensus compounds (n={n})",
                 loc="left", fontsize=13)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig_adme_rules.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_adme_rules.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
