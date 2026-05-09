"""Histogram of the null-model novelty distribution + observed ANPDB marker."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
NM = ROOT / "output" / "null_model"

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def main():
    df = pd.read_csv(NM / "novelty_null_distribution.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, col, observed, label in zip(
        axes,
        ["pass1_rate", "pass2_rate"],
        [11.91, 8.84],
        ["Pass 1 — InChIKey block novel",
         "Pass 2 — block novel + Tanimoto < 0.85"],
    ):
        rate = df[col] * 100
        ax.hist(rate, bins=20, color="#AAAAAA", edgecolor="white", alpha=0.7,
                label=f"Null distribution\n(n={len(df)} bootstrap draws)")
        ax.axvline(observed, color="#D62728", lw=3,
                   label=f"ANPDB observed: {observed:.2f} %")
        ax.text(observed + 0.5, ax.get_ylim()[1] * 0.85,
                f"{observed:.2f} %", color="#D62728", fontsize=12, fontweight="bold")
        ax.axvline(rate.mean(), color="#1C7293", lw=2, ls="--", alpha=0.6,
                   label=f"Null mean: {rate.mean():.2f} %")
        ax.set_xlabel("Novelty rate (%)")
        ax.set_ylabel("Bootstrap draws")
        ax.set_title(label, fontsize=12, loc="left")
        ax.legend(fontsize=10, loc="upper right")

    fig.suptitle(
        "Permutation null model: ANPDB's truly-novel rate is "
        "SIGNIFICANTLY LOWER than expected for a random NP collection",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(NM / "fig_novelty_null.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {NM / 'fig_novelty_null.png'}")


if __name__ == "__main__":
    main()
