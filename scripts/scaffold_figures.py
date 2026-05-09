"""Generate publication figures for Bemis-Murcko scaffold analysis.

Reads precomputed CSVs from output/scaffolds/ and produces:
- fig_scaffold_freq.png        top-20 most frequent scaffolds
- fig_scaffold_diversity.png   cumulative scaffold coverage curve
- fig_generic_freq.png         top-20 generic frameworks
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "scaffolds"

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def main() -> int:
    sc = pd.read_csv(OUT / "scaffold_counts.csv")
    gc = pd.read_csv(OUT / "generic_scaffold_counts.csv")

    # --- Top 20 Murcko scaffolds ---
    top = 20
    df1 = sc.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.35 * top + 1)))
    y = np.arange(len(df1))
    ax.barh(y, df1["n_compounds"], color="#1C7293")
    for yi, (n, frac) in enumerate(zip(df1["n_compounds"], df1["fraction"])):
        ax.text(n + 0.3, yi, f"{100*frac:.1f}%", va="center", fontsize=9, color="#444")
    ax.set_yticks(y)
    labels = df1["scaffold_smiles"].apply(lambda s: s[:55] + "..." if len(s) > 55 else s)
    ax.set_yticklabels(labels, fontsize=8, family="monospace")
    ax.set_xlabel("Number of novel compounds sharing this scaffold")
    ax.set_title("Top 20 Bemis-Murcko scaffolds in the 1,012 novel ANPDB compounds",
                 loc="left", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "fig_scaffold_freq.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_scaffold_freq.png'}")

    # --- Cumulative coverage curve ---
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(1, len(sc) + 1)
    ax.plot(x, sc["cumulative_fraction"], color="#1C7293", lw=2)
    ax.axhline(0.5, color="#D62728", ls="--", lw=1, alpha=0.6)
    n50 = (sc["cumulative_fraction"] >= 0.5).idxmax() + 1
    ax.axvline(n50, color="#D62728", ls="--", lw=1, alpha=0.6)
    ax.text(n50 + 5, 0.52, f"50% coverage at {n50} scaffolds", fontsize=10, color="#D62728")
    ax.set_xlabel("Number of scaffolds (ranked by frequency)")
    ax.set_ylabel("Cumulative fraction of compounds covered")
    ax.set_title("Scaffold diversity: cumulative coverage curve", loc="left", fontsize=13)
    n_total = len(sc)
    n_sing = (sc["n_compounds"] == 1).sum()
    ax.text(0.98, 0.15,
            f"Total unique scaffolds: {n_total:,}\n"
            f"Singletons: {n_sing:,} ({100*n_sing/n_total:.0f}%)\n"
            f"Diversity ratio: {n_total/sc['n_compounds'].sum():.3f}",
            transform=ax.transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#ddd"),
            fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig_scaffold_diversity.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_scaffold_diversity.png'}")

    # --- Top 20 generic frameworks ---
    df2 = gc.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.35 * top + 1)))
    y = np.arange(len(df2))
    ax.barh(y, df2["n_compounds"], color="#2CA02C")
    for yi, (n, frac) in enumerate(zip(df2["n_compounds"], df2["fraction"])):
        ax.text(n + 0.3, yi, f"{100*frac:.1f}%", va="center", fontsize=9, color="#444")
    ax.set_yticks(y)
    labels2 = df2["generic_scaffold_smiles"].apply(lambda s: s[:55] + "..." if len(s) > 55 else s)
    ax.set_yticklabels(labels2, fontsize=8, family="monospace")
    ax.set_xlabel("Number of novel compounds sharing this generic framework")
    ax.set_title("Top 20 generic Murcko frameworks in the 1,012 novel ANPDB compounds",
                 loc="left", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "fig_generic_frameworks.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_generic_frameworks.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
