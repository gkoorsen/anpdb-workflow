"""Deep-dive on the ophthalmology polypharmacology cluster (clusters 4 + 6 from
the IDF-weighted analysis of the 1,012 truly-novel ANPDB compounds).

Outputs (under output/network/ophthalmology/)
---------------------------------------------
- ophthalmology_compounds.csv   per-compound metadata + cluster-level summary
- ophthalmology_targets.csv     per-target compound count + disease bucket
- ophthalmology_top_species.csv plant species frequency in the cluster
- subnetwork_nodes.tsv          Cytoscape-compatible nodes
- subnetwork_edges.tsv          Cytoscape-compatible edges
- figures/
    fig_top_compounds.png       leaderboard by ophthalmology disease reach
    fig_top_targets.png         top hit targets in the cluster
    fig_top_species.png         most-represented source plant species
    fig_disease_buckets.png     buckets (glaucoma, cataract, hearing, other)
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
OUT = ROOT / "output" / "network" / "ophthalmology"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PIDGIN = Path("/Users/gerritkoorsen/PIDGINv4")
HIT_AD60 = PIDGIN / "ANPDB_novel_1012_pair_hits_ad60_with_diseases.tsv"
CLUSTERS = ROOT / "output" / "network" / "compound_clusters_ad60_idf.csv"
ANPDB_RAW = ROOT / "data" / "ANPDB.csv"
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"

OPHTHAL_CLUSTERS = {4, 6}

# Super-category bucketing — collapse near-synonyms
BUCKET_PATTERNS = [
    ("Glaucoma",                 re.compile(r"glaucoma", re.I)),
    ("Cataract / Lens",          re.compile(r"\b(cataract|lens opaci|aphakia|crystalline lens)", re.I)),
    ("Hearing loss / Otology",   re.compile(r"hearing|deaf|otopath|sensorineural", re.I)),
    ("Cornea / Sclera",          re.compile(r"corneal|sclera|hydrophthalm|globe of eye|irido-corneo", re.I)),
    ("Retina / Optic",           re.compile(r"retin|optic|macul|nyctalop", re.I)),
    ("Other eye",                re.compile(r"eye|ocular|ophth", re.I)),
]


def bucket_for(disease: str) -> str | None:
    for label, pat in BUCKET_PATTERNS:
        if pat.search(disease):
            return label
    return None


_sep_re = re.compile(r" {4,}")


def parse_anpdb_metadata() -> pd.DataFrame:
    """Re-parse ANPDB.csv to recover source_species / family / kingdom."""
    rows: list[dict] = []
    with open(ANPDB_RAW, encoding="utf-8") as fh:
        fh.readline()
        for raw in fh:
            fields = _sep_re.split(raw.rstrip("\n"))
            inchi_idx = next(
                (i for i, v in enumerate(fields) if v.startswith("InChI=")), -1
            )
            if inchi_idx < 0:
                continue
            mid = fields[0]
            mol_name = fields[1] if len(fields) > 1 else ""
            # After inchi: source_species, family, kingdom, class_name, subclass_name, synonyms, pmids, links
            tail = fields[inchi_idx + 1:]
            src = tail[0] if len(tail) > 0 else ""
            fam = tail[1] if len(tail) > 1 else ""
            kdm = tail[2] if len(tail) > 2 else ""
            cls = tail[3] if len(tail) > 3 else ""
            sub = tail[4] if len(tail) > 4 else ""
            rows.append({
                "molecule_id": mid, "mol_name": mol_name,
                "source_species": src, "family": fam, "kingdom": kdm,
                "class_name": cls, "subclass_name": sub,
            })
    return pd.DataFrame(rows)


def main() -> int:
    sns.set_theme(context="talk", style="whitegrid", palette="deep")
    pd.set_option("display.max_colwidth", 80)

    print("Loading IDF clusters ...")
    clusters = pd.read_csv(CLUSTERS, dtype={"compound_id": str})
    ophthal_ids = set(clusters.loc[clusters["cluster"].isin(OPHTHAL_CLUSTERS), "compound_id"])
    print(f"  ophthalmology cluster (4+6) compounds: {len(ophthal_ids)}")

    print("Loading hits ...")
    hits = pd.read_csv(HIT_AD60, sep="\t", dtype=str, keep_default_na=False)
    hits["prediction_probability"] = pd.to_numeric(hits["prediction_probability"], errors="coerce")
    hits["disease_set"] = hits["diseases"].apply(
        lambda s: [d.strip() for d in (s or "").split(" | ") if d.strip()]
    )
    sub = hits[hits["compound_id"].isin(ophthal_ids)].copy()
    print(f"  hits in cluster: {len(sub):,}")
    print(f"  cpds: {sub['compound_id'].nunique()}  targets: {sub['uniprot'].nunique()}")

    # bucket per row by mapping each disease in disease_set
    sub["_buckets"] = sub["disease_set"].apply(
        lambda lst: {b for d in lst if (b := bucket_for(d))}
    )
    sub["_any_eye"] = sub["_buckets"].apply(bool)
    eye_hits = sub[sub["_any_eye"]].copy()
    print(f"  eye-relevant hits: {len(eye_hits):,} / {len(sub):,}")

    # --- per-compound summary ---
    print("\nLoading ANPDB metadata + std SMILES ...")
    metadata = parse_anpdb_metadata()
    metadata = metadata[metadata["molecule_id"].isin(ophthal_ids)]
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    smiles_map = dict(zip(novel["molecule_id"], novel["std_smiles"].fillna("")))
    name_map = dict(zip(novel["molecule_id"], novel["mol_name"].fillna("")))

    cmpd_summary: list[dict] = []
    for cid, grp in sub.groupby("compound_id"):
        eye = grp[grp["_any_eye"]]
        bucket_counter = Counter()
        target_set = set()
        max_p_eye = 0.0
        for _, r in eye.iterrows():
            for b in r["_buckets"]:
                bucket_counter[b] += 1
            target_set.add(r["uniprot"])
            max_p_eye = max(max_p_eye, r["prediction_probability"] or 0.0)
        meta = metadata[metadata["molecule_id"] == cid]
        cluster_id = clusters.loc[clusters["compound_id"] == cid, "cluster"].iloc[0]
        cmpd_summary.append({
            "molecule_id": cid,
            "cluster": int(cluster_id),
            "mol_name": meta["mol_name"].iloc[0] if len(meta) else name_map.get(cid, ""),
            "source_species": meta["source_species"].iloc[0] if len(meta) else "",
            "family": meta["family"].iloc[0] if len(meta) else "",
            "class_name": meta["class_name"].iloc[0] if len(meta) else "",
            "subclass_name": meta["subclass_name"].iloc[0] if len(meta) else "",
            "smiles": smiles_map.get(cid, ""),
            "n_eye_targets": len(target_set),
            "n_eye_hits": int(sum(bucket_counter.values())),
            "max_prob_eye": round(max_p_eye, 3),
            "buckets_hit": " | ".join(sorted(bucket_counter)),
            "n_glaucoma_hits": bucket_counter.get("Glaucoma", 0),
            "n_cataract_hits": bucket_counter.get("Cataract / Lens", 0),
            "n_hearing_hits": bucket_counter.get("Hearing loss / Otology", 0),
            "n_cornea_hits": bucket_counter.get("Cornea / Sclera", 0),
        })
    cmpd_df = pd.DataFrame(cmpd_summary).sort_values(
        ["n_eye_targets", "max_prob_eye"], ascending=[False, False]
    )
    cmpd_df.to_csv(OUT / "ophthalmology_compounds.csv", index=False)
    print(f"  wrote {len(cmpd_df)} compound rows")

    # --- top targets ---
    tgt_rows: list[dict] = []
    for tgt, grp in eye_hits.groupby("uniprot"):
        bucket_counter: Counter = Counter()
        for _, r in grp.iterrows():
            for b in r["_buckets"]:
                bucket_counter[b] += 1
        tgt_rows.append({
            "uniprot": tgt,
            "target_label": grp["target_label"].iloc[0],
            "n_compounds": grp["compound_id"].nunique(),
            "n_hits": len(grp),
            "max_prob": round(grp["prediction_probability"].max() or 0.0, 3),
            "median_prob": round(grp["prediction_probability"].median() or 0.0, 3),
            "buckets": " | ".join(f"{b}({n})" for b, n in bucket_counter.most_common()),
        })
    tgt_df = pd.DataFrame(tgt_rows).sort_values(
        ["n_compounds", "max_prob"], ascending=[False, False]
    )
    tgt_df.to_csv(OUT / "ophthalmology_targets.csv", index=False)
    print(f"\nTop 10 ophthalmology-driving targets:")
    print(tgt_df.head(10).to_string(index=False))

    # --- top source species ---
    species_counter: Counter = Counter()
    for src in cmpd_df["source_species"]:
        for s in (src or "").split(" || "):
            s = s.strip()
            if s:
                species_counter[s] += 1
    species_df = pd.DataFrame(
        species_counter.most_common(),
        columns=["source_species", "n_compounds_in_cluster"],
    )
    species_df.to_csv(OUT / "ophthalmology_top_species.csv", index=False)
    print(f"\nTop 10 source species in cluster:")
    print(species_df.head(10).to_string(index=False))

    # --- subnetwork export ---
    nodes: list[dict] = []
    edges: list[dict] = []
    seen = set()

    def add(node_id, ntype, **attrs):
        if node_id in seen:
            return
        seen.add(node_id)
        nodes.append({"id": node_id, "type": ntype, **attrs})

    bucket_meta = {
        "Glaucoma":               "#D62728",
        "Cataract / Lens":        "#FF7F0E",
        "Hearing loss / Otology": "#9467BD",
        "Cornea / Sclera":        "#2CA02C",
        "Retina / Optic":         "#8C564B",
        "Other eye":              "#7F7F7F",
    }
    cmpd_cluster = dict(zip(cmpd_df["molecule_id"], cmpd_df["cluster"]))
    for _, r in eye_hits.iterrows():
        add(r["compound_id"], "compound", cluster=cmpd_cluster.get(r["compound_id"], 0))
        add(r["uniprot"], "target", label=r["target_label"])
        edges.append({
            "source": r["compound_id"], "target": r["uniprot"],
            "type": "compound-target",
            "weight": r["prediction_probability"] or 0.0,
            "ad_percentile": r["ad_percentile"],
        })
        for b in r["_buckets"]:
            add(b, "disease_bucket", color=bucket_meta.get(b, "#444"))
            edges.append({
                "source": r["uniprot"], "target": b,
                "type": "target-disease",
                "weight": 1.0,
            })
    pd.DataFrame(nodes).to_csv(OUT / "subnetwork_nodes.tsv", sep="\t", index=False)
    pd.DataFrame(edges).drop_duplicates().to_csv(
        OUT / "subnetwork_edges.tsv", sep="\t", index=False
    )
    print(f"\nSubnetwork: {len(nodes)} nodes / {len(edges):,} edges")

    # --- figures ---
    # fig: top compounds by # eye targets
    top_n = 20
    df1 = cmpd_df.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.4 * len(df1) + 1)))
    y = np.arange(len(df1))
    ax.barh(y, df1["n_eye_targets"], color="#1C7293")
    for yi, (n_h, mp) in enumerate(zip(df1["n_eye_hits"], df1["max_prob_eye"])):
        ax.text(df1["n_eye_targets"].iloc[yi] + 0.05, yi,
                f"hits={n_h}  max p={mp:.2f}",
                va="center", fontsize=9, color="#444")
    ax.set_yticks(y)
    ax.set_yticklabels(df1["molecule_id"])
    ax.set_xlabel("Distinct eye-relevant targets predicted hit (ad60)")
    ax.set_title(
        "Top compounds in the ophthalmology cluster (n=98)",
        loc="left", fontsize=14,
    )
    fig.tight_layout(); fig.savefig(FIG / "fig_top_compounds.png", dpi=160); plt.close(fig)

    # fig: top targets
    df2 = tgt_df.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.4 * len(df2) + 1)))
    y = np.arange(len(df2))
    ax.barh(y, df2["n_compounds"], color="#EE854A")
    for yi, (label, mp) in enumerate(zip(df2["target_label"], df2["max_prob"])):
        ax.text(df2["n_compounds"].iloc[yi] + 0.5, yi,
                f"{label}  max p={mp:.2f}",
                va="center", fontsize=9, color="#444")
    ax.set_yticks(y)
    ax.set_yticklabels(df2["uniprot"])
    ax.set_xlabel("Distinct ophthalmology-cluster compounds hitting target (ad60)")
    ax.set_title("Targets driving the ophthalmology signal", loc="left", fontsize=14)
    fig.tight_layout(); fig.savefig(FIG / "fig_top_targets.png", dpi=160); plt.close(fig)

    # fig: top species
    df3 = species_df.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, max(5, 0.4 * len(df3) + 1)))
    y = np.arange(len(df3))
    ax.barh(y, df3["n_compounds_in_cluster"], color="#2CA02C")
    ax.set_yticks(y)
    ax.set_yticklabels(df3["source_species"])
    ax.set_xlabel("Distinct ophthalmology-cluster compounds isolated from species")
    ax.set_title("Top source plant species in the ophthalmology cluster",
                 loc="left", fontsize=14)
    fig.tight_layout(); fig.savefig(FIG / "fig_top_species.png", dpi=160); plt.close(fig)

    # fig: disease buckets
    bucket_totals: Counter = Counter()
    bucket_cmpds: dict[str, set[str]] = {}
    for _, r in eye_hits.iterrows():
        for b in r["_buckets"]:
            bucket_totals[b] += 1
            bucket_cmpds.setdefault(b, set()).add(r["compound_id"])
    bucket_df = pd.DataFrame([
        {"bucket": b, "n_hits": bucket_totals[b], "n_compounds": len(bucket_cmpds[b])}
        for b in bucket_totals
    ]).sort_values("n_hits", ascending=False)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(bucket_df) + 1)))
    y = np.arange(len(bucket_df))
    ax.barh(y, bucket_df["n_hits"], color=[bucket_meta.get(b, "#444") for b in bucket_df["bucket"]])
    ax.set_yticks(y)
    ax.set_yticklabels(bucket_df["bucket"])
    for yi, n_c in enumerate(bucket_df["n_compounds"]):
        ax.text(bucket_df["n_hits"].iloc[yi] + 0.5, yi,
                f"{n_c} cpds", va="center", fontsize=10, color="#444")
    ax.invert_yaxis()
    ax.set_xlabel("Hit count by ophthalmology disease bucket")
    ax.set_title("Ophthalmology disease bucket breakdown", loc="left", fontsize=14)
    fig.tight_layout(); fig.savefig(FIG / "fig_disease_buckets.png", dpi=160); plt.close(fig)

    print(f"\nAll outputs in {OUT}/")
    print("\nDisease bucket totals:")
    print(bucket_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
