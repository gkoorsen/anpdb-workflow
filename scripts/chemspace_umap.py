"""Stage 2: UMAP dimensionality reduction + chemical space plot.

Reads precomputed .npy fingerprint arrays from output/chemspace/
and produces a UMAP scatter plot showing novel ANPDB compounds
vs a COCONUT background sample.

Output: output/chemspace/fig_chemical_space_umap.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "chemspace"


def main() -> int:
    print("Loading fingerprint arrays ...")
    novel_fps = np.load(OUT / "novel_fps.npy")
    coconut_fps = np.load(OUT / "coconut_sample_fps.npy")
    labels = pd.read_csv(OUT / "labels.csv")
    print(f"  novel: {novel_fps.shape}  coconut: {coconut_fps.shape}")

    X = np.vstack([novel_fps, coconut_fps]).astype(np.float32)
    sources = labels["source"].values

    print("Running UMAP (n_neighbors=30, min_dist=0.3, metric=jaccard) ...")
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric="jaccard",
        n_components=2,
        random_state=42,
        n_jobs=1,
    )
    embedding = reducer.fit_transform(X)
    print(f"  embedding shape: {embedding.shape}")

    coconut_mask = sources == "COCONUT"
    novel_mask = sources == "Novel ANPDB"

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(
        embedding[coconut_mask, 0], embedding[coconut_mask, 1],
        s=4, alpha=0.15, c="#AAAAAA", label=f"COCONUT sample (n={coconut_mask.sum():,})",
        rasterized=True,
    )
    ax.scatter(
        embedding[novel_mask, 0], embedding[novel_mask, 1],
        s=8, alpha=0.5, c="#D62728", label=f"Novel ANPDB (n={novel_mask.sum():,})",
        rasterized=True,
    )
    ax.legend(fontsize=11, markerscale=3)
    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(
        "Chemical space: 1,012 novel ANPDB compounds vs COCONUT",
        fontsize=14, loc="left",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "fig_chemical_space_umap.png", dpi=200)
    plt.close(fig)
    print(f"Saved figure to {OUT / 'fig_chemical_space_umap.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
