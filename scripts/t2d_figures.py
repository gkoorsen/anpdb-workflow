"""Visual summary of ANPDB compound -> target -> Type 2 diabetes associations.

Input compound pool: ANPDB compounds that passed Lipinski's rule of five, PAINS,
and ADMET pre-filters before PIDGINv4 target prediction. Every compound shown
is already drug-likeness-viable; the figures focus purely on T2D relevance.

Reads the PIDGINv4 hits files (ad60 = looser, ad90 = stringent) and produces:
  fig1_t2d_targets.png       targets implicated in T2D, ranked by compound hits
  fig2_compound_target_heatmap.png   prediction-probability heatmap of top T2D
                                     hits (compounds x targets) at ad60
  fig3_proba_vs_ad.png       prediction_probability vs applicability domain
                             scatter, colored by target, faceted by AD cutoff
  fig4_top_compounds.png     leaderboard of compounds with the most T2D-target
                             hits at ad60, annotated with their ad90 promotion
                             status
  fig5_t2d_disease_breadth.png   compounds vs distinct T2D-related diseases
                                 they hit (ad60)
  t2d_hits_ad60.csv / t2d_hits_ad90.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PIDGIN = Path("/Users/gerritkoorsen/PIDGINv4")
HITS60 = PIDGIN / "FILE_ADMET_SMALL_FIXED_FROM_ANPDB_small_pair_hits_ad60_with_diseases.tsv"
HITS90 = PIDGIN / "FILE_ADMET_SMALL_FIXED_FROM_ANPDB_small_pair_hits_ad90_with_diseases.tsv"

T2D_STRENGTH_MIN = 4  # legacy strength cutoff (no longer the primary filter)
USE_STRICT_FILTER = True  # option B: keep only targets with NIDDM term named
T2D_STRICT = "Diabetes Mellitus, Non-Insulin-Dependent"
T2D_RELATED = {
    "Diabetes Mellitus, Non-Insulin-Dependent",
    "Insulin Resistance",
    "Hyperglycemia",
    "Hyperinsulinism",
    "Diabetic Nephropathy",
    "Diabetic Neuropathies",
    "Diabetic Retinopathy",
    "Diabetic Angiopathies",
    "Complications of Diabetes Mellitus",
    "MICROVASCULAR COMPLICATIONS OF DIABETES, SUSCEPTIBILITY TO, 1(finding)",
    "MICROVASCULAR COMPLICATIONS OF DIABETES, SUSCEPTIBILITY TO, 3 (finding)",
    "Gestational Diabetes",
}

sns.set_theme(context="talk", style="whitegrid", palette="deep")


def load_hits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df["prediction_probability"] = pd.to_numeric(df["prediction_probability"], errors="coerce")
    df["ad_percentile"] = pd.to_numeric(df["ad_percentile"], errors="coerce")
    df["disease_set"] = df["diseases"].apply(
        lambda s: {d.strip() for d in s.split(" | ") if d.strip()} if s else set()
    )
    df["t2d_related"] = df["disease_set"].apply(lambda s: bool(s & T2D_RELATED))
    df["t2d_strict"] = df["disease_set"].apply(lambda s: T2D_STRICT in s)
    df["n_t2d_diseases_on_target"] = df["disease_set"].apply(lambda s: len(s & T2D_RELATED))
    df["n_total_diseases_on_target"] = df["disease_set"].apply(len)
    # T2D-association strength: count of T2D-related diseases on the target,
    # plus a +2 bonus when the strict NIDDM term is named explicitly.
    df["t2d_strength"] = (
        df["n_t2d_diseases_on_target"] + 2 * df["t2d_strict"].astype(int)
    )
    df["t2d_share"] = df["n_t2d_diseases_on_target"] / df["n_total_diseases_on_target"].clip(lower=1)
    return df


def target_strength_table(df: pd.DataFrame) -> pd.DataFrame:
    """One row per target with T2D strength metrics + compound counts."""
    if df.empty:
        return pd.DataFrame()
    g = (
        df.groupby("uniprot")
        .agg(
            n_compounds=("compound_id", "nunique"),
            t2d_strength=("t2d_strength", "first"),
            t2d_share=("t2d_share", "first"),
            n_t2d_diseases=("n_t2d_diseases_on_target", "first"),
            n_total_diseases=("n_total_diseases_on_target", "first"),
            t2d_strict=("t2d_strict", "first"),
            median_prob=("prediction_probability", "median"),
        )
        .reset_index()
    )
    return g.sort_values(["t2d_strength", "n_compounds"], ascending=[False, False])


def fig1_targets(t60: pd.DataFrame, t90: pd.DataFrame, novel_ids: set[str], path: Path) -> None:
    """Bar chart: targets ranked by number of distinct compound hits."""
    def per_target(df: pd.DataFrame) -> pd.DataFrame:
        g = (
            df.groupby("uniprot")
            .agg(
                n_hits=("compound_id", "size"),
                n_compounds=("compound_id", "nunique"),
                n_novel_compounds=(
                    "compound_id",
                    lambda s: s[s.isin(novel_ids)].nunique(),
                ),
                median_prob=("prediction_probability", "median"),
            )
            .reset_index()
        )
        return g.sort_values("n_compounds", ascending=False)

    g60 = per_target(t60)
    g90 = per_target(t90)

    if g60.empty:
        return

    n_show = min(15, len(g60))
    top = g60.head(n_show).iloc[::-1]
    promoted_ad90 = set(g90["uniprot"]) if not g90.empty else set()

    fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * n_show + 1)))
    ypos = np.arange(len(top))
    ax.barh(ypos, top["n_compounds"], color="#4878D0", label="distinct compounds (ad60)")
    ax.barh(
        ypos,
        top["n_novel_compounds"],
        color="#EE854A",
        label="of which in tightened novel set",
    )
    for y, u in zip(ypos, top["uniprot"]):
        if u in promoted_ad90:
            ax.text(
                top.set_index("uniprot").loc[u, "n_compounds"] + 0.5,
                y,
                "* ad90",
                va="center",
                fontsize=10,
                color="#444",
            )
    ax.set_yticks(ypos)
    ax.set_yticklabels(top["uniprot"])
    ax.set_xlabel("Distinct ANPDB compounds predicted to hit (ad60)")
    ax.set_title(
        "T2D-related targets — compounds pre-filtered (Lipinski + PAINS + ADMET)\n"
        f"({len(g60)} targets at ad60, {len(g90)} at ad90 — * marks promotion to ad90)",
        loc="left",
        fontsize=14,
    )
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig2_heatmap(t60: pd.DataFrame, novel_ids: set[str], path: Path) -> None:
    """Heatmap of prediction probability for top compounds x top targets."""
    if t60.empty:
        return

    top_targets = (
        t60.groupby("uniprot")["compound_id"].nunique().sort_values(ascending=False).head(12).index
    )
    df = t60[t60["uniprot"].isin(top_targets)].copy()
    cmpd_score = (
        df.groupby("compound_id")["prediction_probability"].max().sort_values(ascending=False)
    )
    top_cmpds = cmpd_score.head(30).index
    df = df[df["compound_id"].isin(top_cmpds)]
    pivot = df.pivot_table(
        index="compound_id",
        columns="uniprot",
        values="prediction_probability",
        aggfunc="max",
    )
    pivot = pivot.reindex(index=top_cmpds, columns=top_targets)

    annot_idx = [
        (cid + "  ●") if cid in novel_ids else cid for cid in pivot.index
    ]

    fig, ax = plt.subplots(figsize=(10, max(6, 0.32 * len(pivot) + 2)))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="viridis",
        vmin=0.5,
        vmax=1.0,
        cbar_kws={"label": "max prediction probability"},
        linewidths=0.3,
        linecolor="#fff",
    )
    ax.set_yticklabels(annot_idx, rotation=0, fontsize=10)
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_xlabel("UniProt target")
    ax.set_ylabel("ANPDB compound (● = in 1,370 novel set)")
    ax.set_title(
        "Top ANPDB compound × T2D-related target prediction probabilities (ad60)",
        loc="left",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig3_proba_vs_ad(t60: pd.DataFrame, t90: pd.DataFrame, path: Path) -> None:
    """Scatter: prediction_probability vs ad_percentile for T2D-target hits."""
    if t60.empty:
        return
    df = t60.copy()
    df["passes_ad90"] = df.set_index(["compound_id", "uniprot"]).index.isin(
        t90.set_index(["compound_id", "uniprot"]).index
    )
    fig, ax = plt.subplots(figsize=(9, 6.5))
    sns.scatterplot(
        data=df,
        x="ad_percentile",
        y="prediction_probability",
        hue="passes_ad90",
        palette={True: "#D65F5F", False: "#5C6B73"},
        s=70,
        alpha=0.75,
        edgecolor="white",
        ax=ax,
    )
    ax.axvline(60, ls="--", lw=1, color="#999", label="ad60 cutoff")
    ax.axvline(90, ls="--", lw=1, color="#D65F5F", label="ad90 cutoff")
    ax.axhline(0.5, ls=":", lw=1, color="#999")
    ax.set_xlim(55, 101)
    ax.set_ylim(0.45, 1.02)
    ax.set_xlabel("Applicability-domain percentile")
    ax.set_ylabel("Prediction probability")
    ax.set_title(
        f"All T2D-related hits at ad60 (n={len(df):,}) — red = also passes ad90",
        loc="left",
        fontsize=13,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig4_top_compounds(t60: pd.DataFrame, t90: pd.DataFrame, novel_ids: set[str], path: Path) -> None:
    """Leaderboard: compounds ranked by # distinct T2D targets they hit at ad60."""
    if t60.empty:
        return
    g = (
        t60.groupby("compound_id")
        .agg(
            n_targets=("uniprot", "nunique"),
            n_t2d_strict=("t2d_strict", "sum"),
            max_prob=("prediction_probability", "max"),
            median_prob=("prediction_probability", "median"),
        )
        .reset_index()
    )
    promoted = set(zip(t90["compound_id"], t90["uniprot"]))
    g["n_targets_ad90"] = g["compound_id"].apply(
        lambda c: len(set(t for cc, t in promoted if cc == c))
    )
    g["is_novel"] = g["compound_id"].isin(novel_ids)
    g = g.sort_values(["n_targets", "max_prob"], ascending=[False, False]).head(20)

    fig, ax = plt.subplots(figsize=(11, max(5, 0.45 * len(g) + 1)))
    ypos = np.arange(len(g))[::-1]
    ax.barh(
        ypos,
        g["n_targets"],
        color=["#EE854A" if n else "#4878D0" for n in g["is_novel"]],
        label="targets at ad60",
    )
    ax.barh(
        ypos,
        g["n_targets_ad90"],
        color="#3F4D63",
        label="of which also at ad90",
    )
    labels = [
        f"{cid}{'  ●' if novel else ''}"
        for cid, novel in zip(g["compound_id"], g["is_novel"])
    ]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    for y, (cid, mp) in zip(ypos, zip(g["compound_id"], g["max_prob"])):
        ax.text(
            g.set_index("compound_id").loc[cid, "n_targets"] + 0.05,
            y,
            f"max p={mp:.2f}",
            va="center",
            fontsize=9,
            color="#444",
        )
    ax.set_xlabel("Distinct T2D-related targets predicted hit")
    ax.set_title(
        "Top ANPDB compounds by T2D-target breadth (orange = in 1,370 novel set)",
        loc="left",
        fontsize=13,
    )
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def fig5_disease_breadth(t60: pd.DataFrame, novel_ids: set[str], path: Path) -> None:
    """For each compound, count distinct T2D-related diseases its targets touch."""
    if t60.empty:
        return
    rows: list[dict] = []
    for cid, sub in t60.groupby("compound_id"):
        diseases: set[str] = set()
        for s in sub["disease_set"]:
            diseases |= s & T2D_RELATED
        rows.append(
            {
                "compound_id": cid,
                "t2d_diseases_touched": len(diseases),
                "n_targets": sub["uniprot"].nunique(),
                "is_novel": cid in novel_ids,
            }
        )
    g = pd.DataFrame(rows)
    g = g[g["t2d_diseases_touched"] > 0].sort_values(
        ["t2d_diseases_touched", "n_targets"], ascending=[False, False]
    ).head(20)

    fig, ax = plt.subplots(figsize=(10, max(5, 0.42 * len(g) + 1)))
    ypos = np.arange(len(g))[::-1]
    ax.barh(
        ypos,
        g["t2d_diseases_touched"],
        color=["#EE854A" if n else "#4878D0" for n in g["is_novel"]],
    )
    labels = [
        f"{cid}{'  ●' if novel else ''}"
        for cid, novel in zip(g["compound_id"], g["is_novel"])
    ]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Distinct T2D-related diseases reachable via predicted targets (ad60)")
    ax.set_title("Compounds with broadest T2D-related disease reach", loc="left", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> int:
    print(f"Loading hits ...", file=sys.stderr)
    h60 = load_hits(HITS60)
    h90 = load_hits(HITS90)

    t60 = h60[h60["t2d_related"]].copy()
    t90 = h90[h90["t2d_related"]].copy()
    print(
        f"  T2D-related hits at ad60: {len(t60):,} pairs / "
        f"{t60['compound_id'].nunique():,} compounds / {t60['uniprot'].nunique():,} targets",
        file=sys.stderr,
    )
    print(
        f"  T2D-related hits at ad90: {len(t90):,} pairs / "
        f"{t90['compound_id'].nunique():,} compounds / {t90['uniprot'].nunique():,} targets",
        file=sys.stderr,
    )
    strict60 = h60[h60["t2d_strict"]]
    strict90 = h90[h90["t2d_strict"]]
    print(
        f"  STRICT 'Diabetes Mellitus, Non-Insulin-Dependent' hits at ad60: "
        f"{len(strict60):,} pairs / {strict60['compound_id'].nunique():,} compounds",
        file=sys.stderr,
    )
    print(
        f"  STRICT T2D hits at ad90: {len(strict90):,} pairs / "
        f"{strict90['compound_id'].nunique():,} compounds",
        file=sys.stderr,
    )

    # Use the TIGHTENED novel set: standardization + Tanimoto < 0.85 vs COCONUT.
    novel_path = ROOT / "output" / "truly_novel_molecule_ids_std.txt"
    novel_ids = set(open(novel_path).read().split())
    print(f"  novel-set size: {len(novel_ids):,} (from {novel_path.name})", file=sys.stderr)

    # Per-target T2D-association strength (informational; not the primary filter)
    target_table = target_strength_table(t60)
    target_table.to_csv(ROOT / "output" / "t2d_target_strength.csv", index=False)
    print("\nPer-target T2D strength (top 10):", file=sys.stderr)
    print(target_table.head(10).to_string(index=False), file=sys.stderr)

    # Primary filter: option B — keep targets with strict NIDDM term (t2d_strict=True)
    if USE_STRICT_FILTER:
        kept = set(target_table.loc[target_table["t2d_strict"], "uniprot"])
        filter_label = "t2d_strict=True (option B)"
    else:
        kept = set(target_table.loc[target_table["t2d_strength"] >= T2D_STRENGTH_MIN, "uniprot"])
        filter_label = f"t2d_strength >= {T2D_STRENGTH_MIN}"
    t60s = t60[t60["uniprot"].isin(kept)].copy()
    t90s = t90[t90["uniprot"].isin(kept)].copy()
    print(
        f"\nAfter filter ({filter_label}):"
        f"\n  targets retained: {len(kept)} of {len(target_table)}"
        f"\n  ad60 hits: {len(t60s):,} pairs / {t60s['compound_id'].nunique():,} compounds"
        f"\n    of which novel: {t60s['compound_id'].isin(novel_ids).sum()} pairs / "
        f"{t60s.loc[t60s['compound_id'].isin(novel_ids),'compound_id'].nunique()} compounds"
        f"\n  ad90 hits: {len(t90s):,} pairs / {t90s['compound_id'].nunique():,} compounds",
        file=sys.stderr,
    )

    t60.drop(columns=["disease_set"]).to_csv(ROOT / "output" / "t2d_hits_ad60.csv", index=False)
    t90.drop(columns=["disease_set"]).to_csv(ROOT / "output" / "t2d_hits_ad90.csv", index=False)
    t60s.drop(columns=["disease_set"]).to_csv(ROOT / "output" / "t2d_hits_ad60_strong.csv", index=False)
    t90s.drop(columns=["disease_set"]).to_csv(ROOT / "output" / "t2d_hits_ad90_strong.csv", index=False)

    # Broad-T2D figures (unchanged baseline)
    fig1_targets(t60, t90, novel_ids, OUT / "fig1_t2d_targets.png")
    fig2_heatmap(t60, novel_ids, OUT / "fig2_compound_target_heatmap.png")
    fig3_proba_vs_ad(t60, t90, OUT / "fig3_proba_vs_ad.png")
    fig4_top_compounds(t60, t90, novel_ids, OUT / "fig4_top_compounds.png")
    fig5_disease_breadth(t60, novel_ids, OUT / "fig5_t2d_disease_breadth.png")

    # Strength-filtered figures
    fig0_target_strength(target_table, OUT / "fig0_t2d_target_strength.png")
    fig1_targets(t60s, t90s, novel_ids, OUT / "fig1s_t2d_targets_strong.png")
    fig2_heatmap(t60s, novel_ids, OUT / "fig2s_compound_target_heatmap_strong.png")
    fig3_proba_vs_ad(t60s, t90s, OUT / "fig3s_proba_vs_ad_strong.png")
    fig4_top_compounds(t60s, t90s, novel_ids, OUT / "fig4s_top_compounds_strong.png")
    fig5_disease_breadth(t60s, novel_ids, OUT / "fig5s_t2d_disease_breadth_strong.png")

    print(f"\nFigures written to {OUT}/", file=sys.stderr)
    return 0


def fig0_target_strength(table: pd.DataFrame, path: Path) -> None:
    """Per-target T2D-association strength: bars colored by share, marked
    with strict-NIDDM badge, sized by compound hit count.
    """
    if table.empty:
        return
    df = table.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, max(4, 0.45 * len(df) + 1)))
    ypos = np.arange(len(df))
    cmap = plt.get_cmap("YlOrRd")
    colors = [cmap(min(0.15 + s * 6, 0.95)) for s in df["t2d_share"]]
    ax.barh(ypos, df["t2d_strength"], color=colors, edgecolor="#333", linewidth=0.4)

    labels = [
        f"{u}{' *' if strict else ''}"
        for u, strict in zip(df["uniprot"], df["t2d_strict"])
    ]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)

    for y, (s, c, share) in enumerate(zip(df["t2d_strength"], df["n_compounds"], df["t2d_share"])):
        ax.text(
            s + 0.08,
            y,
            f"{c} cpds   share={share*100:.1f}%",
            va="center",
            fontsize=9,
            color="#333",
        )

    ax.axvline(T2D_STRENGTH_MIN, color="#D65F5F", ls="--", lw=1.2,
               label=f"strength filter (>= {T2D_STRENGTH_MIN})")
    ax.set_xlabel("T2D-association strength (T2D-related disease count + 2 if strict NIDDM term)")
    ax.set_title(
        "Per-target T2D-association strength — '*' = strict NIDDM term in target's diseases\n"
        "bar shade = share of target's diseases that are T2D-related (darker = more specific)",
        loc="left",
        fontsize=12,
    )
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
