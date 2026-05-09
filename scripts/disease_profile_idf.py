"""IDF-weighted disease-profile analysis for the 1,012 novel ANPDB compounds.

Rationale
---------
Disease terms inherited from PIDGIN target-disease annotations are uneven:
- Promiscuous terms (e.g. "Autosomal recessive predisposition", "Schizophrenia")
  appear on many target disease lists and dominate the raw leaderboard.
- Specific terms appear on a handful of targets and are more informative.

We treat each target as a "document" and each disease as a "term", and
re-weight (compound, disease) scores by the disease's IDF:

    IDF(d) = log(N_targets / df(d))

where df(d) = number of distinct targets with disease d in their annotation,
and N_targets is the size of the PIDGIN target panel (the disease lookup we
built from FILE_ADMET_SMALL_FIXED... pair_metrics).

Per (compound, disease) IDF-weighted score:

    s(c, d) = max_t in T_hit(c) [ prob(c,t) * I(d in diseases(t)) ] * IDF(d)

where T_hit(c) = set of targets compound c hits at the chosen AD/probability cutoffs.

Outputs (under output/network/)
-------------------------------
- disease_idf.csv                        IDF score per disease
- disease_leaderboard_ad60_idf.csv       leaderboard re-ranked by IDF reach
- compound_disease_matrix_ad60_idf.csv   IDF-weighted matrix used for clustering
- compound_clusters_ad60_idf.csv         IDF-clustered compound assignments
- figures/fig_idf_distribution.png       per-disease IDF distribution (sanity)
- figures/fig_disease_leaderboard_ad60_idf.png
- figures/fig_cluster_summary_ad60_idf.png
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "network"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PIDGIN = Path("/Users/gerritkoorsen/PIDGINv4")
HIT_AD60 = PIDGIN / "ANPDB_novel_1012_pair_hits_ad60_with_diseases.tsv"
LOOKUP = PIDGIN / "disease_lookup_target_level.csv"

K_CLUSTERS = 6
TOP_DISEASES = 30

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def compute_idf(lookup_path: Path) -> tuple[pd.Series, int]:
    """IDF over the full PIDGIN target-disease lookup."""
    csv.field_size_limit(1 << 30)
    targets_with: dict[str, set[str]] = {}
    n_targets = 0
    with open(lookup_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            n_targets += 1
            tgt_id = (row["uniprot"], row["activity_threshold"])
            for d in (row["diseases"] or "").split(" | "):
                d = d.strip()
                if not d:
                    continue
                targets_with.setdefault(d, set()).add(tgt_id)
    df = pd.Series(
        {d: len(t) for d, t in targets_with.items()},
        name="df",
    )
    idf = pd.Series(
        {d: math.log(n_targets / dfv) for d, dfv in df.items()},
        name="idf",
    )
    out = pd.concat([df, idf], axis=1).sort_values("idf", ascending=False)
    out.index.name = "disease"
    return out, n_targets


def load_hits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df["prediction_probability"] = pd.to_numeric(df["prediction_probability"], errors="coerce")
    df["disease_set"] = df["diseases"].apply(
        lambda s: [d.strip() for d in (s or "").split(" | ") if d.strip()]
    )
    return df


def explode(hits: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in hits.iterrows():
        for d in r["disease_set"]:
            rows.append({
                "compound_id": r["compound_id"],
                "uniprot": r["uniprot"],
                "disease": d,
                "prob": r["prediction_probability"],
            })
    return pd.DataFrame(rows)


def fig_idf(idf_table: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(idf_table["idf"], bins=40, color="#1C7293", edgecolor="white")
    ax.set_xlabel("IDF (higher = more specific)")
    ax.set_ylabel("Number of disease terms")
    ax.set_title("Per-disease IDF over the PIDGIN target panel", loc="left", fontsize=13)
    ax.text(
        0.98, 0.95,
        f"Total disease terms : {len(idf_table):,}\n"
        f"Median IDF          : {idf_table['idf'].median():.2f}\n"
        f"Promiscuous terms (IDF<1): {(idf_table['idf']<1).sum():,}",
        transform=ax.transAxes, ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="#ddd"),
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_leaderboard_idf(lb: pd.DataFrame, path: Path, top: int = 25) -> None:
    df = lb.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.4 * len(df) + 1)))
    y = np.arange(len(df))
    ax.barh(y, df["weighted_reach"], color="#1C7293")
    for yi, (n_c, idf) in enumerate(zip(df["n_compounds"], df["idf"])):
        ax.text(
            df["weighted_reach"].iloc[yi] * 1.02 + 0.05,
            yi,
            f"cpds={n_c}  idf={idf:.2f}",
            va="center",
            fontsize=9,
            color="#444",
        )
    ax.set_yticks(y)
    ax.set_yticklabels([d if len(d) <= 60 else d[:57] + "..." for d in df["disease"]])
    ax.set_xlabel("Weighted reach = sum over novel compounds of (max prob x IDF)")
    ax.set_title(
        f"Top {top} diseases by IDF-weighted reach (ad60 hits)",
        loc="left", fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig_cluster_summary(matrix: pd.DataFrame, clusters: pd.DataFrame, path: Path) -> None:
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
        ax.set_title(f"Cluster {cid}  (n={n_cpds} compounds) — mean IDF-weighted score per disease",
                     loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> int:
    print("Computing IDF over PIDGIN target panel ...")
    idf_table, n_targets = compute_idf(LOOKUP)
    print(f"  N targets={n_targets:,}  N disease terms={len(idf_table):,}")
    print(f"  median IDF={idf_table['idf'].median():.2f}")
    print(f"  IDF<1 (promiscuous): {(idf_table['idf']<1).sum():,}")
    idf_table.to_csv(OUT / "disease_idf.csv")
    fig_idf(idf_table, FIG / "fig_idf_distribution.png")

    print("\nLoading ad60 hits ...")
    hits = load_hits(HIT_AD60)
    long = explode(hits)
    print(f"  hits={len(hits):,}  triples={len(long):,}  cpds={hits['compound_id'].nunique():,}")

    # join IDF
    long = long.merge(
        idf_table["idf"].reset_index(),
        left_on="disease", right_on="disease", how="left",
    )
    long["idf"] = long["idf"].fillna(idf_table["idf"].median())
    long["weighted"] = long["prob"] * long["idf"]

    # Per (compound, disease) score = max weighted across all hit targets
    cd = (
        long.groupby(["compound_id", "disease"])
        .agg(weighted=("weighted", "max"),
             max_prob=("prob", "max"),
             idf=("idf", "first"))
        .reset_index()
    )

    # Disease leaderboard (IDF-weighted)
    lb = (
        cd.groupby("disease")
        .agg(
            n_compounds=("compound_id", "nunique"),
            weighted_reach=("weighted", "sum"),
            mean_weighted=("weighted", "mean"),
            max_prob=("max_prob", "max"),
            idf=("idf", "first"),
        )
        .reset_index()
        .sort_values("weighted_reach", ascending=False)
    )
    lb.to_csv(OUT / "disease_leaderboard_ad60_idf.csv", index=False)
    fig_leaderboard_idf(lb, FIG / "fig_disease_leaderboard_ad60_idf.png")
    print(f"\nTop 15 diseases by IDF-weighted reach:")
    for _, r in lb.head(15).iterrows():
        print(f"  {r['disease'][:60]:<60}  cpds={r['n_compounds']:>3}  "
              f"idf={r['idf']:.2f}  reach={r['weighted_reach']:.2f}  max p={r['max_prob']:.2f}")

    # IDF-weighted compound x top-disease matrix
    top_diseases = lb.head(TOP_DISEASES)["disease"].tolist()
    matrix = (
        cd[cd["disease"].isin(top_diseases)]
        .pivot_table(index="compound_id", columns="disease", values="weighted", aggfunc="max")
        .reindex(columns=top_diseases)
        .fillna(0.0)
    )
    matrix.to_csv(OUT / "compound_disease_matrix_ad60_idf.csv")
    print(f"\nMatrix shape (cpds x top-{TOP_DISEASES}-diseases): {matrix.shape}")

    # Re-cluster on IDF-weighted profiles
    if len(matrix) >= 2:
        Z = linkage(pdist(matrix.values, metric="cosine"), method="average")
        cluster_ids = fcluster(Z, t=K_CLUSTERS, criterion="maxclust")
        clusters_df = pd.DataFrame({"compound_id": matrix.index, "cluster": cluster_ids})
        clusters_df.to_csv(OUT / "compound_clusters_ad60_idf.csv", index=False)
        fig_cluster_summary(matrix, clusters_df, FIG / "fig_cluster_summary_ad60_idf.png")
        print(f"\nCluster sizes (IDF-weighted, k={K_CLUSTERS}):")
        for cid, n in clusters_df["cluster"].value_counts().sort_index().items():
            print(f"  cluster {cid}: {n} compounds")

    print(f"\nAll outputs in {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
