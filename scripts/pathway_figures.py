"""Generate pathway enrichment figure from g:Profiler results.

Reads output/pathways/enrichment_results.csv and produces:
- fig_pathway_enrichment.png    dot plot of enriched terms by source
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "pathways"

sns.set_theme(context="talk", style="whitegrid", palette="deep")

SOURCE_COLORS = {"KEGG": "#1C7293", "REAC": "#EE854A", "GO:BP": "#2CA02C"}


def main() -> int:
    df = pd.read_csv(OUT / "enrichment_results.csv")
    df["-log10(p)"] = -np.log10(df["p_value"].clip(lower=1e-20))
    df = df.sort_values("-log10(p)", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(df) + 1)))
    y = np.arange(len(df))
    colors = [SOURCE_COLORS.get(s, "#444") for s in df["source"]]
    sizes = df["intersection_size"] * 40

    ax.scatter(df["-log10(p)"], y, c=colors, s=sizes, alpha=0.7, edgecolors="white", lw=0.5)
    ax.set_yticks(y)
    term_labels = [
        f"[{row['source']}] {row['term_name'][:50]}" + ("..." if len(row['term_name']) > 50 else "")
        for _, row in df.iterrows()
    ]
    ax.set_yticklabels(term_labels, fontsize=9)

    for yi, (_, row) in enumerate(df.iterrows()):
        ax.text(row["-log10(p)"] + 0.05, yi,
                f"{row['intersection_size']}/{row['term_size']}",
                va="center", fontsize=8, color="#444")

    ax.axvline(-np.log10(0.05), color="#D62728", ls="--", lw=1, alpha=0.5)
    ax.text(-np.log10(0.05) + 0.02, len(df) - 0.5, "p=0.05", fontsize=9, color="#D62728")

    ax.set_xlabel("-log10(p-value)", fontsize=12)
    ax.set_title("Pathway enrichment of 21 consensus targets\n(g:Profiler, g:SCS correction)",
                 loc="left", fontsize=13)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=10, label=s)
        for s, c in SOURCE_COLORS.items()
    ]
    ax.legend(handles=legend_elements, fontsize=10, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "fig_pathway_enrichment.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {OUT / 'fig_pathway_enrichment.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
