"""Disease-profile analysis for the 1,012 truly-novel ANPDB compounds.

Inputs
------
- /Users/gerritkoorsen/PIDGINv4/ANPDB_novel_1012_pair_hits_ad60_with_diseases.tsv
- /Users/gerritkoorsen/PIDGINv4/ANPDB_novel_1012_pair_hits_ad90_with_diseases.tsv
- /Users/gerritkoorsen/anpdb-novelty/output/anpdb_truly_novel_std.csv   (for SMILES)

Outputs (under output/network/)
-------------------------------
- disease_leaderboard_ad60.csv          ranked diseases with reach metrics
- disease_leaderboard_ad90.csv
- compound_disease_matrix_ad60.csv      compound x top-disease matrix (max prob)
- compound_disease_matrix_ad90.csv
- compound_clusters_ad60.csv            cluster id + top diseases per compound
- network_nodes.tsv / network_edges.tsv Cytoscape-compatible compound-target-disease graph (ad60)
- figures:
    fig_disease_leaderboard.png
    fig_compound_disease_heatmap.png
    fig_cluster_dendrogram.png
    fig_cluster_summary.png

Conventions
-----------
- "hit" = is_hit==1 AND meets_ad_cutoff==1
- Disease terms exploded from the per-target "diseases" list ("|"-separated).
- Per-(compound, disease) score = max prediction_probability over all hit
  targets that mention that disease.
- IDF weighting penalizes diseases that appear on many targets (less specific).
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "network"
OUT.mkdir(parents=True, exist_ok=True)
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)

PIDGIN = Path("/Users/gerritkoorsen/PIDGINv4")
HIT_AD60 = PIDGIN / "ANPDB_novel_1012_pair_hits_ad60_with_diseases.tsv"
HIT_AD90 = PIDGIN / "ANPDB_novel_1012_pair_hits_ad90_with_diseases.tsv"

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def load_hits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df["prediction_probability"] = pd.to_numeric(df["prediction_probability"], errors="coerce")
    df["ad_percentile"] = pd.to_numeric(df["ad_percentile"], errors="coerce")
    df["disease_set"] = df["diseases"].apply(
        lambda s: [d.strip() for d in s.split(" | ") if d.strip()] if s else []
    )
    return df


def explode_disease_pairs(hits: pd.DataFrame) -> pd.DataFrame:
    """Long-form: one row per (compound, target, disease) tuple."""
    rows: list[dict] = []
    for _, r in hits.iterrows():
        for d in r["disease_set"]:
            rows.append(
                {
                    "compound_id": r["compound_id"],
                    "uniprot": r["uniprot"],
                    "disease": d,
                    "prediction_probability": r["prediction_probability"],
                    "ad_percentile": r["ad_percentile"],
                }
            )
    return pd.DataFrame(rows)


def disease_leaderboard(long: pd.DataFrame) -> pd.DataFrame:
    g = (
        long.groupby("disease")
        .agg(
            n_hits=("compound_id", "size"),
            n_compounds=("compound_id", "nunique"),
            n_targets=("uniprot", "nunique"),
            mean_prob=("prediction_probability", "mean"),
            max_prob=("prediction_probability", "max"),
        )
        .reset_index()
        .sort_values(["n_compounds", "n_hits"], ascending=[False, False])
    )
    return g


def compound_disease_matrix(long: pd.DataFrame, top_diseases: list[str]) -> pd.DataFrame:
    """Compound x disease pivot with max probability per cell."""
    sub = long[long["disease"].isin(top_diseases)]
    pivot = sub.pivot_table(
        index="compound_id",
        columns="disease",
        values="prediction_probability",
        aggfunc="max",
    ).fillna(0.0)
    pivot = pivot.reindex(columns=top_diseases).fillna(0.0)
    return pivot


def fig_leaderboard(lb: pd.DataFrame, label: str, path: Path, top: int = 25) -> None:
    df = lb.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.4 * len(df) + 1)))
    y = np.arange(len(df))
    ax.barh(y, df["n_compounds"], color="#1C7293", label="distinct novel compounds")
    for yi, (n_c, n_h, mp) in enumerate(
        zip(df["n_compounds"], df["n_hits"], df["max_prob"])
    ):
        ax.text(
            n_c + 0.5,
            yi,
            f"hits={n_h}  max p={mp:.2f}",
            va="center",
            fontsize=9,
            color="#444",
        )
    ax.set_yticks(y)
    ax.set_yticklabels([d if len(d) <= 60 else d[:57] + "..." for d in df["disease"]])
    ax.set_xlabel("Distinct novel compounds reaching the disease (via predicted hit targets)")
    ax.set_title(
        f"Top {top} diseases reachable from the 1,012 truly-novel ANPDB compounds ({label})",
        loc="left",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_heatmap(matrix: pd.DataFrame, label: str, path: Path,
                top_compounds: int = 40, top_diseases: int = 25) -> None:
    cmpd_score = matrix.sum(axis=1).sort_values(ascending=False)
    rows = cmpd_score.head(top_compounds).index
    cols = matrix.columns[:top_diseases]
    sub = matrix.loc[rows, cols]
    fig, ax = plt.subplots(figsize=(max(10, 0.32 * len(cols) + 4),
                                    max(6, 0.28 * len(rows) + 1)))
    sns.heatmap(
        sub,
        ax=ax,
        cmap="viridis",
        vmin=0.5,
        vmax=1.0,
        cbar_kws={"label": "max prediction probability"},
        linewidths=0.3,
        linecolor="#fff",
    )
    ax.set_xticks(np.arange(len(sub.columns)) + 0.5)
    ax.set_xticklabels(
        [d if len(d) <= 35 else d[:33] + "..." for d in sub.columns],
        rotation=45, ha="right", fontsize=9,
    )
    ax.set_yticks(np.arange(len(sub.index)) + 0.5)
    ax.set_yticklabels(sub.index, rotation=0, fontsize=9)
    ax.set_title(
        f"Top compound x top disease prediction probabilities ({label})",
        loc="left",
        fontsize=13,
    )
    ax.set_xlabel("")
    ax.set_ylabel("ANPDB novel compound")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def cluster_compounds(matrix: pd.DataFrame, k: int = 6) -> pd.DataFrame:
    if len(matrix) < 2:
        return pd.DataFrame()
    X = matrix.values
    # Cosine distance on disease profile vectors
    Z = linkage(pdist(X, metric="cosine"), method="average")
    cluster_ids = fcluster(Z, t=k, criterion="maxclust")
    out = pd.DataFrame(
        {"compound_id": matrix.index, "cluster": cluster_ids}
    )
    return out, Z


def fig_dendrogram(matrix: pd.DataFrame, Z, k: int, path: Path) -> None:
    if len(matrix) < 2:
        return
    fig, ax = plt.subplots(figsize=(14, max(6, 0.18 * len(matrix) + 1)))
    dendrogram(
        Z,
        labels=matrix.index.tolist(),
        leaf_rotation=0,
        leaf_font_size=7,
        orientation="left",
        color_threshold=0,
        ax=ax,
    )
    ax.set_xlabel("cosine distance (compound disease profiles)")
    ax.set_title(
        f"Compound clustering by predicted-disease profile (k={k} clusters)",
        loc="left",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_cluster_summary(matrix: pd.DataFrame, clusters: pd.DataFrame, path: Path) -> None:
    """For each cluster, show its top 8 dominant diseases as a stacked bar chart."""
    df = matrix.copy()
    df["cluster"] = df.index.map(dict(zip(clusters["compound_id"], clusters["cluster"])))
    summary: dict[int, pd.Series] = {}
    for cid, sub in df.groupby("cluster"):
        col_means = sub.drop(columns="cluster").mean()
        summary[int(cid)] = col_means.sort_values(ascending=False).head(8)
    n_clusters = len(summary)
    fig, axes = plt.subplots(n_clusters, 1, figsize=(11, 2.2 * n_clusters), sharex=False)
    if n_clusters == 1:
        axes = [axes]
    for ax, (cid, top) in zip(axes, sorted(summary.items())):
        n_cpds = (clusters["cluster"] == cid).sum()
        ax.barh(np.arange(len(top)), top.values[::-1], color="#1C7293")
        ax.set_yticks(np.arange(len(top)))
        ax.set_yticklabels(
            [d if len(d) <= 55 else d[:53] + "..." for d in top.index[::-1]],
            fontsize=9,
        )
        ax.set_xlim(0, 1.05)
        ax.set_title(f"Cluster {cid}  (n={n_cpds} compounds) — mean prob per disease",
                     loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def export_cytoscape(hits: pd.DataFrame, lb: pd.DataFrame,
                     long: pd.DataFrame, top_diseases: list[str]) -> None:
    """Write nodes.tsv + edges.tsv suitable for Cytoscape."""
    nodes: list[dict] = []
    edges: list[dict] = []
    seen_node = set()

    def add(node_id: str, ntype: str, **attrs):
        if node_id in seen_node:
            return
        seen_node.add(node_id)
        nodes.append({"id": node_id, "type": ntype, **attrs})

    # restrict to top diseases for tractability
    sub = long[long["disease"].isin(top_diseases)]

    for _, r in hits.iterrows():
        add(r["compound_id"], "compound")
        add(r["uniprot"], "target", target_label=r.get("target_label", ""))
        edges.append({
            "source": r["compound_id"],
            "target": r["uniprot"],
            "type": "compound-target",
            "weight": r["prediction_probability"],
            "ad_percentile": r["ad_percentile"],
        })
    # one target-disease edge per (target, disease) pair
    seen_td: set[tuple[str, str]] = set()
    for _, r in sub.iterrows():
        td = (r["uniprot"], r["disease"])
        if td in seen_td:
            continue
        seen_td.add(td)
        add(r["disease"], "disease")
        edges.append({
            "source": r["uniprot"],
            "target": r["disease"],
            "type": "target-disease",
            "weight": 1.0,
        })

    pd.DataFrame(nodes).to_csv(OUT / "network_nodes.tsv", sep="\t", index=False)
    pd.DataFrame(edges).to_csv(OUT / "network_edges.tsv", sep="\t", index=False)


def run(label: str, path: Path) -> None:
    print(f"\n=== {label} ===")
    hits = load_hits(path)
    print(f"  hits     : {len(hits):,}")
    print(f"  cpds     : {hits['compound_id'].nunique():,}")
    print(f"  targets  : {hits['uniprot'].nunique():,}")

    long = explode_disease_pairs(hits)
    print(f"  (cpd,target,disease) tuples: {len(long):,}")

    lb = disease_leaderboard(long)
    lb.to_csv(OUT / f"disease_leaderboard_{label}.csv", index=False)
    print(f"  diseases : {len(lb):,}")
    print(f"  top 10 diseases by compound coverage:")
    for _, r in lb.head(10).iterrows():
        print(f"    {r['disease'][:60]:<60}  cpds={r['n_compounds']:>3}  hits={r['n_hits']:>3}  max p={r['max_prob']:.2f}")

    # Top-30 diseases for matrix and figures
    top_diseases = lb.head(30)["disease"].tolist()
    matrix = compound_disease_matrix(long, top_diseases)
    matrix.to_csv(OUT / f"compound_disease_matrix_{label}.csv")
    print(f"  matrix shape: {matrix.shape}")

    fig_leaderboard(lb, label, FIG / f"fig_disease_leaderboard_{label}.png")
    fig_heatmap(matrix, label, FIG / f"fig_compound_disease_heatmap_{label}.png")

    if label == "ad60":
        # cluster only on ad60 (more breadth)
        clusters_df, Z = cluster_compounds(matrix, k=6)
        clusters_df.to_csv(OUT / "compound_clusters_ad60.csv", index=False)
        fig_dendrogram(matrix, Z, k=6, path=FIG / "fig_cluster_dendrogram_ad60.png")
        fig_cluster_summary(matrix, clusters_df, FIG / "fig_cluster_summary_ad60.png")
        # Cytoscape export from ad60
        export_cytoscape(hits, lb, long, top_diseases)
        print(f"  cytoscape: nodes={len(set(hits['compound_id'])|set(hits['uniprot'])|set(top_diseases)):,}")


def main() -> int:
    run("ad60", HIT_AD60)
    run("ad90", HIT_AD90)
    print(f"\nAll outputs in {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
