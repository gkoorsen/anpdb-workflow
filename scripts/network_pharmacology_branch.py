"""Disease-target / PPI / hub-gene branch for ANPDB lead hypotheses.

This script adds the conventional network-pharmacology branch that is missing
from the compound-first PIDGIN + ChEMBL-NN workflow:

1. fetch disease-associated targets from Open Targets by disease ID;
2. combine the disease module with configured ANPDB consensus drug targets;
3. fetch a human STRING PPI network;
4. rank network hubs and compute drug-target proximity to the disease module;
5. optionally run g:Profiler enrichment on overlap/near-neighbour module genes.

The branch is deliberately disease-specific. It is meant for hypotheses such as
CYP1B1/glaucoma, SGLT2/T2D, and MAO-B/Parkinson disease, not for one large
all-diseases network.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
import time
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from urllib import error, parse, request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "network_pharmacology" / "lead_hypotheses.csv"
TARGET_MAP = ROOT / "configs" / "network_pharmacology" / "target_map.csv"
CONSENSUS = ROOT / "output" / "chembl_nn" / "consensus_pidgin_chembl_nn_ad60.tsv"
OUT = ROOT / "output" / "network_pharmacology"

OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
STRING_URL = "https://string-db.org/api/tsv/network"
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

CALLER_ID = "anpdb-novelty-network-pharmacology"
RANDOM_SEED = 42

MPL_CACHE = OUT / ".matplotlib-cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(context="talk", style="whitegrid", palette="deep")


@dataclass
class Hypothesis:
    hypothesis: str
    disease_name: str
    disease_id: str
    lead_target_symbol: str
    lead_target_uniprot: str
    top_disease_targets: int
    min_open_targets_score: float
    string_required_score: int
    string_add_nodes: int
    random_iterations: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run disease-target/PPI/hub-gene network pharmacology branch."
    )
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--target-map", type=Path, default=TARGET_MAP)
    parser.add_argument("--consensus", type=Path, default=CONSENSUS)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument(
        "--hypothesis",
        action="append",
        help="Run only this hypothesis key. Repeat to run several.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Refetch API data even when cached output files already exist.",
    )
    parser.add_argument(
        "--skip-enrichment",
        action="store_true",
        help="Skip g:Profiler enrichment on the final module genes.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and print planned jobs without contacting APIs.",
    )
    return parser.parse_args()


def read_hypotheses(path: Path) -> list[Hypothesis]:
    rows = pd.read_csv(path, dtype=str, keep_default_na=False)
    hypotheses: list[Hypothesis] = []
    for _, row in rows.iterrows():
        hypotheses.append(
            Hypothesis(
                hypothesis=row["hypothesis"],
                disease_name=row["disease_name"],
                disease_id=row.get("disease_id") or row.get("efo_id"),
                lead_target_symbol=row["lead_target_symbol"],
                lead_target_uniprot=row["lead_target_uniprot"],
                top_disease_targets=int(row.get("top_disease_targets") or 300),
                min_open_targets_score=float(row.get("min_open_targets_score") or 0.05),
                string_required_score=int(row.get("string_required_score") or 700),
                string_add_nodes=int(row.get("string_add_nodes") or 0),
                random_iterations=int(row.get("random_iterations") or 1000),
            )
        )
    return hypotheses


def read_target_map(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df = df[df.get("include_in_network", "1").astype(str) != "0"].copy()
    return {
        row["uniprot"]: {"symbol": row["symbol"], "name": row.get("name", "")}
        for _, row in df.iterrows()
        if row.get("uniprot") and row.get("symbol")
    }


def post_json(url: str, payload: dict, timeout: int = 90, tries: int = 3) -> dict:
    data = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json", "User-Agent": CALLER_ID}
    last_exc: Exception | None = None
    for attempt in range(1, tries + 1):
        try:
            req = request.Request(url, data=data, headers=headers, method="POST")
            with request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_exc = exc
            if attempt < tries:
                time.sleep(2 * attempt)
    raise RuntimeError(f"POST failed for {url}: {last_exc}") from last_exc


def search_open_targets_disease(query_string: str) -> tuple[str, str] | None:
    query = """
    query diseaseSearch($queryString: String!) {
      search(
        queryString: $queryString,
        entityNames: ["disease"],
        page: {index: 0, size: 10}
      ) {
        hits {
          id
          entity
          name
          description
        }
      }
    }
    """
    payload = {"query": query, "variables": {"queryString": query_string}}
    response = post_json(OPEN_TARGETS_URL, payload)
    if response.get("errors"):
        raise RuntimeError(json.dumps(response["errors"], indent=2))
    hits = response.get("data", {}).get("search", {}).get("hits") or []
    disease_hits = [hit for hit in hits if hit.get("id") and hit.get("name")]
    if not disease_hits:
        return None
    hit = disease_hits[0]
    return str(hit["id"]), str(hit["name"])


def fetch_open_targets_for_id(disease_id: str, size: int) -> tuple[pd.DataFrame, str] | None:
    query = """
    query diseaseTargets($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: 0, size: $size}) {
          count
          rows {
            score
            target {
              id
              approvedSymbol
              approvedName
              proteinIds {
                id
                source
              }
            }
          }
        }
      }
    }
    """
    payload = {
        "query": query,
        "variables": {"efoId": disease_id, "size": size},
    }
    response = post_json(OPEN_TARGETS_URL, payload)
    if response.get("errors"):
        raise RuntimeError(json.dumps(response["errors"], indent=2))
    disease = response.get("data", {}).get("disease")
    if not disease:
        return None

    rows = []
    for item in disease["associatedTargets"]["rows"]:
        target = item["target"]
        uniprots = [
            p["id"]
            for p in target.get("proteinIds") or []
            if p.get("source") == "uniprot_swissprot"
        ]
        rows.append(
            {
                "disease_id": disease["id"],
                "disease_name": disease["name"],
                "target_ensembl": target.get("id", ""),
                "symbol": target.get("approvedSymbol", ""),
                "target_name": target.get("approvedName", ""),
                "uniprot": "|".join(uniprots),
                "open_targets_score": item.get("score", 0.0),
            }
        )
    df = pd.DataFrame(rows)
    if len(df):
        df["open_targets_score"] = pd.to_numeric(
            df["open_targets_score"], errors="coerce"
        ).fillna(0.0)
        df = df.sort_values("open_targets_score", ascending=False)
    return df, str(disease["name"])


def fetch_open_targets(h: Hypothesis) -> tuple[pd.DataFrame, str]:
    fetched = fetch_open_targets_for_id(h.disease_id, h.top_disease_targets)
    if fetched is None:
        resolved = search_open_targets_disease(h.disease_name)
        if resolved is None:
            raise RuntimeError(
                f"Open Targets returned no disease for {h.disease_id} "
                f"and search found no match for {h.disease_name!r}"
            )
        resolved_id, resolved_name = resolved
        print(
            f"  resolved {h.disease_name!r} via Open Targets search: "
            f"{resolved_name} ({resolved_id})"
        )
        fetched = fetch_open_targets_for_id(resolved_id, h.top_disease_targets)
        if fetched is None:
            raise RuntimeError(
                f"Open Targets search resolved {h.disease_name!r} to {resolved_id}, "
                "but disease target lookup still failed."
            )
    df, disease_name = fetched
    if len(df):
        df = df[df["open_targets_score"] >= h.min_open_targets_score].copy()
    return df, disease_name


def cached_open_targets(h: Hypothesis, out_dir: Path, force: bool) -> tuple[pd.DataFrame, str]:
    path = out_dir / "disease_targets_opentargets.csv"
    if path.exists() and not force:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "open_targets_score" in df:
            df["open_targets_score"] = pd.to_numeric(
                df["open_targets_score"], errors="coerce"
            ).fillna(0.0)
        disease_name = df["disease_name"].iloc[0] if len(df) else h.disease_name
        return df, disease_name

    df, disease_name = fetch_open_targets(h)
    df.to_csv(path, index=False)
    return df, disease_name


def build_drug_targets(
    h: Hypothesis,
    consensus_path: Path,
    target_map: dict[str, dict[str, str]],
) -> pd.DataFrame:
    base = {
        "symbol": h.lead_target_symbol,
        "uniprot": h.lead_target_uniprot,
        "name": target_map.get(h.lead_target_uniprot, {}).get("name", ""),
        "role": "lead_target",
        "supporting_compounds": 0,
        "supporting_compound_ids": "",
    }
    if not consensus_path.exists():
        return pd.DataFrame([base])

    consensus = pd.read_csv(consensus_path, sep="\t", dtype=str, keep_default_na=False)
    lead_compounds = sorted(
        consensus.loc[
            consensus["uniprot"] == h.lead_target_uniprot, "compound_id"
        ].unique()
    )
    base["supporting_compounds"] = len(lead_compounds)
    base["supporting_compound_ids"] = "|".join(lead_compounds)

    rows = [base]
    if lead_compounds:
        co_hits = consensus[consensus["compound_id"].isin(lead_compounds)].copy()
        for uniprot, sub in co_hits.groupby("uniprot"):
            if uniprot == h.lead_target_uniprot or uniprot not in target_map:
                continue
            rows.append(
                {
                    "symbol": target_map[uniprot]["symbol"],
                    "uniprot": uniprot,
                    "name": target_map[uniprot].get("name", ""),
                    "role": "mapped_consensus_cotarget",
                    "supporting_compounds": sub["compound_id"].nunique(),
                    "supporting_compound_ids": "|".join(sorted(sub["compound_id"].unique())),
                }
            )
    return pd.DataFrame(rows).drop_duplicates("symbol")


def fetch_string_edges(symbols: list[str], h: Hypothesis) -> pd.DataFrame:
    if not symbols:
        return pd.DataFrame()
    data = parse.urlencode(
        {
            "identifiers": "\r".join(sorted(set(symbols))),
            "species": "9606",
            "required_score": str(h.string_required_score),
            "network_type": "functional",
            "add_nodes": str(h.string_add_nodes),
            "caller_identity": CALLER_ID,
        }
    ).encode("utf-8")
    headers = {"User-Agent": CALLER_ID}
    req = request.Request(STRING_URL, data=data, headers=headers, method="POST")
    with request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    if not text.strip():
        return pd.DataFrame()
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    return pd.DataFrame(rows)


def normalize_string_edges(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame(columns=["source", "target", "score"])
    source_col = "preferredName_A" if "preferredName_A" in raw.columns else "stringId_A"
    target_col = "preferredName_B" if "preferredName_B" in raw.columns else "stringId_B"
    score_col = "score" if "score" in raw.columns else "combined_score"
    out = raw[[source_col, target_col, score_col]].copy()
    out.columns = ["source", "target", "score"]
    out["score"] = pd.to_numeric(out["score"], errors="coerce").fillna(0.0)
    if out["score"].max() > 1.0:
        out["score"] = out["score"] / 1000.0
    out = out[out["source"] != out["target"]].drop_duplicates(["source", "target"])
    return out.sort_values("score", ascending=False)


def cached_string_edges(
    symbols: list[str],
    h: Hypothesis,
    out_dir: Path,
    force: bool,
) -> pd.DataFrame:
    path = out_dir / "ppi_edges_string.tsv"
    if path.exists() and not force:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    raw = fetch_string_edges(symbols, h)
    edges = normalize_string_edges(raw)
    edges.to_csv(path, sep="\t", index=False)
    return edges


def build_adjacency(edges: pd.DataFrame) -> dict[str, dict[str, float]]:
    adj: dict[str, dict[str, float]] = defaultdict(dict)
    for _, row in edges.iterrows():
        a = str(row["source"])
        b = str(row["target"])
        score = float(row.get("score", 0.0) or 0.0)
        if not a or not b or a == b:
            continue
        adj[a][b] = max(score, adj[a].get(b, 0.0))
        adj[b][a] = max(score, adj[b].get(a, 0.0))
    return adj


def bfs_distances(adj: dict[str, dict[str, float]], sources: set[str]) -> dict[str, int]:
    dist: dict[str, int] = {}
    q: deque[str] = deque()
    for s in sources:
        if s in adj and s not in dist:
            dist[s] = 0
            q.append(s)
    while q:
        v = q.popleft()
        for w in adj.get(v, {}):
            if w not in dist:
                dist[w] = dist[v] + 1
                q.append(w)
    return dist


def betweenness_centrality(adj: dict[str, dict[str, float]]) -> dict[str, float]:
    nodes = list(adj)
    cb = {v: 0.0 for v in nodes}
    for s in nodes:
        stack: list[str] = []
        pred = {w: [] for w in nodes}
        sigma = dict.fromkeys(nodes, 0.0)
        sigma[s] = 1.0
        dist = dict.fromkeys(nodes, -1)
        dist[s] = 0
        q: deque[str] = deque([s])
        while q:
            v = q.popleft()
            stack.append(v)
            for w in adj[v]:
                if dist[w] < 0:
                    q.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = dict.fromkeys(nodes, 0.0)
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w]:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                cb[w] += delta[w]
    for v in cb:
        cb[v] /= 2.0
    return cb


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce").fillna(0.0)
    std = values.std(ddof=0)
    if std == 0 or math.isnan(std):
        return values * 0.0
    return (values - values.mean()) / std


def graph_metrics(
    edges: pd.DataFrame,
    disease_targets: pd.DataFrame,
    drug_targets: pd.DataFrame,
    lead_symbol: str,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    adj = build_adjacency(edges)
    disease_symbols = set(disease_targets["symbol"].dropna().astype(str))
    drug_symbols = set(drug_targets["symbol"].dropna().astype(str))
    all_symbols = set(adj) | disease_symbols | drug_symbols

    ot_score = {
        str(r["symbol"]): float(r.get("open_targets_score", 0.0) or 0.0)
        for _, r in disease_targets.iterrows()
    }
    degree = {n: len(adj.get(n, {})) for n in all_symbols}
    weighted_degree = {n: sum(adj.get(n, {}).values()) for n in all_symbols}
    betweenness = betweenness_centrality(adj) if len(adj) <= 1200 else {}
    dist_to_disease = bfs_distances(adj, disease_symbols)
    dist_to_drug = bfs_distances(adj, drug_symbols)

    rows = []
    for symbol in sorted(all_symbols):
        is_disease = symbol in disease_symbols
        is_drug = symbol in drug_symbols
        role = "ppi_neighbor"
        if is_disease and is_drug:
            role = "overlap_drug_disease_target"
        elif is_drug:
            role = "predicted_drug_target"
        elif is_disease:
            role = "disease_target"
        if symbol == lead_symbol:
            role = "lead_target" if not is_disease else "lead_target_disease_overlap"
        rows.append(
            {
                "symbol": symbol,
                "role": role,
                "is_disease_target": int(is_disease),
                "is_predicted_drug_target": int(is_drug),
                "is_lead_target": int(symbol == lead_symbol),
                "open_targets_score": ot_score.get(symbol, 0.0),
                "degree": degree.get(symbol, 0),
                "weighted_degree": weighted_degree.get(symbol, 0.0),
                "betweenness": betweenness.get(symbol, 0.0),
                "distance_to_nearest_disease_target": dist_to_disease.get(symbol, np.nan),
                "distance_to_nearest_drug_target": dist_to_drug.get(symbol, np.nan),
            }
        )
    df = pd.DataFrame(rows)
    df["hub_score"] = (
        zscore(df["weighted_degree"])
        + zscore(df["betweenness"])
        + zscore(df["open_targets_score"])
    )
    return df.sort_values(["hub_score", "weighted_degree"], ascending=False), adj


def proximity_test(
    adj: dict[str, dict[str, float]],
    disease_symbols: set[str],
    drug_symbols: set[str],
    iterations: int,
) -> dict[str, float | int | str]:
    nodes = set(adj)
    disease_present = disease_symbols & nodes
    drug_present = drug_symbols & nodes
    dist_to_disease = bfs_distances(adj, disease_present)
    observed_distances = [
        dist_to_disease[s] for s in drug_present if s in dist_to_disease
    ]
    observed = float(np.mean(observed_distances)) if observed_distances else np.nan

    result: dict[str, float | int | str] = {
        "n_network_nodes": len(nodes),
        "n_network_edges": int(sum(len(v) for v in adj.values()) / 2),
        "n_disease_targets": len(disease_symbols),
        "n_disease_targets_in_ppi": len(disease_present),
        "n_drug_targets": len(drug_symbols),
        "n_drug_targets_in_ppi": len(drug_present),
        "observed_mean_min_distance": observed,
        "random_iterations": iterations,
        "empirical_p_leq_observed": np.nan,
        "drug_targets_in_ppi": "|".join(sorted(drug_present)),
    }
    if not drug_present or not disease_present or math.isnan(observed):
        return result

    candidates = sorted(nodes - disease_present)
    k = len(drug_present)
    if len(candidates) < k:
        candidates = sorted(nodes)
    rng = random.Random(RANDOM_SEED)
    random_means = []
    for _ in range(iterations):
        sample = set(rng.sample(candidates, k))
        sample_distances = [dist_to_disease[s] for s in sample if s in dist_to_disease]
        if sample_distances:
            random_means.append(float(np.mean(sample_distances)))
    if random_means:
        leq = sum(1 for value in random_means if value <= observed)
        result["empirical_p_leq_observed"] = (leq + 1) / (len(random_means) + 1)
        result["random_mean_min_distance_mean"] = float(np.mean(random_means))
        result["random_mean_min_distance_sd"] = float(np.std(random_means))
    return result


def run_gprofiler(genes: list[str]) -> pd.DataFrame:
    if not genes:
        return pd.DataFrame()
    payload = {
        "organism": "hsapiens",
        "query": sorted(set(genes)),
        "sources": ["KEGG", "REAC", "GO:BP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }
    response = post_json(GPROFILER_URL, payload)
    rows = []
    for r in response.get("result", []):
        rows.append(
            {
                "source": r.get("source", ""),
                "term_id": r.get("native", ""),
                "term_name": r.get("name", ""),
                "p_value": r.get("p_value", 1.0),
                "term_size": r.get("term_size", 0),
                "intersection_size": r.get("intersection_size", 0),
                "query_size": r.get("query_size", 0),
                "intersections": " | ".join(map(str, r.get("intersections", []))),
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["source", "p_value"])


def module_genes_for_enrichment(nodes: pd.DataFrame, limit: int = 100) -> list[str]:
    sub = nodes[
        (nodes["is_predicted_drug_target"] == 1)
        | (nodes["role"] == "overlap_drug_disease_target")
        | (nodes["distance_to_nearest_drug_target"].fillna(99) <= 1)
    ].copy()
    if len(sub) < 5:
        sub = nodes.head(limit).copy()
    return sub.head(limit)["symbol"].dropna().astype(str).tolist()


def write_summary(
    h: Hypothesis,
    disease_name: str,
    disease_id: str,
    drug_targets: pd.DataFrame,
    disease_targets: pd.DataFrame,
    nodes: pd.DataFrame,
    proximity: dict[str, float | int | str],
    out_dir: Path,
    enrichment: pd.DataFrame,
) -> None:
    overlap = nodes[nodes["role"].isin(["overlap_drug_disease_target", "lead_target_disease_overlap"])]
    top_hubs = nodes.head(10)
    lead = nodes[nodes["symbol"] == h.lead_target_symbol]
    lead_ot = float(lead["open_targets_score"].iloc[0]) if len(lead) else 0.0
    lead_dist = lead["distance_to_nearest_disease_target"].iloc[0] if len(lead) else np.nan
    pval = proximity.get("empirical_p_leq_observed", np.nan)

    lines = [
        f"# {h.hypothesis}",
        "",
        f"Disease module: **{disease_name}** (`{disease_id}`)",
        f"Lead target: **{h.lead_target_symbol}** (`{h.lead_target_uniprot}`)",
        "",
        "## Counts",
        "",
        f"- Open Targets disease targets retained: {len(disease_targets):,}",
        f"- Drug targets in this branch: {len(drug_targets):,}",
        f"- STRING nodes / edges: {proximity.get('n_network_nodes', 0):,} / {proximity.get('n_network_edges', 0):,}",
        f"- Drug targets present in STRING network: {proximity.get('drug_targets_in_ppi', '') or 'none'}",
        "",
        "## Disease-Network Position",
        "",
        f"- Direct drug/disease overlap targets: {', '.join(overlap['symbol'].tolist()) if len(overlap) else 'none'}",
        f"- Lead Open Targets score: {lead_ot:.3f}",
        f"- Lead shortest PPI distance to disease module: {lead_dist if not pd.isna(lead_dist) else 'not connected'}",
        f"- Mean minimum distance, drug targets to disease module: {proximity.get('observed_mean_min_distance')}",
        f"- Empirical randomization p-value (lower/closer is better): {pval}",
        "",
        "## Top Hub/Context Genes",
        "",
        "| rank | symbol | role | Open Targets score | degree | weighted degree | betweenness |",
        "|---:|---|---|---:|---:|---:|---:|",
    ]
    for rank, (_, row) in enumerate(top_hubs.iterrows(), start=1):
        lines.append(
            f"| {rank} | {row['symbol']} | {row['role']} | "
            f"{float(row['open_targets_score']):.3f} | {int(row['degree'])} | "
            f"{float(row['weighted_degree']):.2f} | {float(row['betweenness']):.2f} |"
        )
    if len(enrichment):
        lines += [
            "",
            "## Top Enrichment Terms",
            "",
            "| source | term | p-value | hits |",
            "|---|---|---:|---:|",
        ]
        for _, row in enrichment.head(10).iterrows():
            lines.append(
                f"| {row['source']} | {str(row['term_name'])[:80]} | "
                f"{float(row['p_value']):.2e} | {int(row['intersection_size'])}/{int(row['term_size'])} |"
            )
    lines += [
        "",
        "## Interpretation Guardrail",
        "",
        "This branch provides disease-network context. It does not prove that the ANPDB compound acts through the hub genes unless those genes are also predicted or experimentally validated compound targets.",
        "",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def fig_hubs(nodes: pd.DataFrame, out_dir: Path, title: str) -> None:
    top = nodes.head(20).iloc[::-1].copy()
    palette = {
        "lead_target": "#D62728",
        "lead_target_disease_overlap": "#D62728",
        "overlap_drug_disease_target": "#9467BD",
        "predicted_drug_target": "#FF7F0E",
        "disease_target": "#1C7293",
        "ppi_neighbor": "#7F7F7F",
    }
    colors = [palette.get(r, "#7F7F7F") for r in top["role"]]
    fig, ax = plt.subplots(figsize=(10, max(5, 0.35 * len(top) + 1)))
    ax.barh(np.arange(len(top)), top["hub_score"], color=colors)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["symbol"], fontsize=9)
    ax.set_xlabel("Composite hub score")
    ax.set_title(title, loc="left", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_hub_candidates.png", dpi=160)
    plt.close(fig)


def write_index(results: list[dict], out_dir: Path) -> None:
    if not results:
        return
    df = pd.DataFrame(results)
    df.to_csv(out_dir / "summary_index.csv", index=False)
    lines = [
        "# Network Pharmacology Summary Index",
        "",
        "| hypothesis | disease | lead target | direct overlap | lead OT score | STRING nodes/edges | proximity p | top context genes |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for _, row in df.iterrows():
        direct_overlap = str(row["direct_overlap"]).replace("|", ", ")
        top_context = str(row["top_context_genes"]).replace("|", ", ")
        lines.append(
            f"| {row['hypothesis']} | {row['disease_name']} | {row['lead_target']} | "
            f"{direct_overlap} | {float(row['lead_open_targets_score']):.3f} | "
            f"{int(row['string_nodes'])}/{int(row['string_edges'])} | "
            f"{float(row['empirical_p_leq_observed']):.3g} | {top_context} |"
        )
    (out_dir / "summary_index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_one(
    h: Hypothesis,
    args: argparse.Namespace,
    target_map: dict[str, dict[str, str]],
) -> dict:
    out_dir = args.out / h.hypothesis
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {h.hypothesis} ===")
    print(f"  disease: {h.disease_name} ({h.disease_id})")
    print(f"  lead target: {h.lead_target_symbol} ({h.lead_target_uniprot})")

    drug_targets = build_drug_targets(h, args.consensus, target_map)
    drug_targets.to_csv(out_dir / "drug_targets.csv", index=False)
    print(f"  drug targets: {len(drug_targets)}")

    if args.dry_run:
        print("  dry-run: skipping Open Targets, STRING and g:Profiler")
        return {"hypothesis": h.hypothesis, "dry_run": True}

    disease_targets, disease_name = cached_open_targets(h, out_dir, args.force)
    resolved_disease_id = (
        str(disease_targets["disease_id"].iloc[0]) if len(disease_targets) else h.disease_id
    )
    print(f"  disease targets retained: {len(disease_targets)}")

    disease_symbols = set(disease_targets["symbol"].dropna().astype(str))
    drug_symbols = set(drug_targets["symbol"].dropna().astype(str))
    string_inputs = sorted(disease_symbols | drug_symbols)
    edges = cached_string_edges(string_inputs, h, out_dir, args.force)
    print(f"  STRING edges: {len(edges)}")

    nodes, adj = graph_metrics(edges, disease_targets, drug_targets, h.lead_target_symbol)
    nodes.to_csv(out_dir / "ppi_nodes_ranked.csv", index=False)
    hub_candidates = nodes.head(50)
    hub_candidates.to_csv(out_dir / "hub_gene_candidates.csv", index=False)

    overlap = nodes[nodes["role"].isin(["overlap_drug_disease_target", "lead_target_disease_overlap"])]
    overlap.to_csv(out_dir / "overlap_targets.csv", index=False)

    proximity = proximity_test(adj, disease_symbols, drug_symbols, h.random_iterations)
    pd.DataFrame([proximity]).to_csv(out_dir / "network_proximity.csv", index=False)
    print(
        "  proximity:",
        f"observed={proximity.get('observed_mean_min_distance')}",
        f"p={proximity.get('empirical_p_leq_observed')}",
    )

    enrichment = pd.DataFrame()
    if not args.skip_enrichment:
        genes = module_genes_for_enrichment(nodes)
        enrichment = run_gprofiler(genes)
        enrichment.to_csv(out_dir / "enrichment_results.csv", index=False)
        print(f"  enrichment terms: {len(enrichment)}")

    fig_hubs(nodes, out_dir, f"{h.hypothesis}: STRING hub/context genes")
    write_summary(
        h,
        disease_name,
        resolved_disease_id,
        drug_targets,
        disease_targets,
        nodes,
        proximity,
        out_dir,
        enrichment,
    )
    print(f"  wrote {out_dir}")
    lead = nodes[nodes["symbol"] == h.lead_target_symbol]
    overlap_symbols = overlap["symbol"].dropna().astype(str).tolist()
    return {
        "hypothesis": h.hypothesis,
        "disease_name": disease_name,
        "disease_id": resolved_disease_id,
        "lead_target": h.lead_target_symbol,
        "lead_uniprot": h.lead_target_uniprot,
        "direct_overlap": "|".join(overlap_symbols) if overlap_symbols else "",
        "lead_open_targets_score": float(lead["open_targets_score"].iloc[0]) if len(lead) else 0.0,
        "string_nodes": proximity.get("n_network_nodes", 0),
        "string_edges": proximity.get("n_network_edges", 0),
        "observed_mean_min_distance": proximity.get("observed_mean_min_distance", np.nan),
        "empirical_p_leq_observed": proximity.get("empirical_p_leq_observed", np.nan),
        "top_context_genes": "|".join(nodes.head(5)["symbol"].astype(str).tolist()),
        "enrichment_terms": len(enrichment),
        "out_dir": str(out_dir),
    }


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    hypotheses = read_hypotheses(args.config)
    if args.hypothesis:
        wanted = set(args.hypothesis)
        hypotheses = [h for h in hypotheses if h.hypothesis in wanted]
        missing = wanted - {h.hypothesis for h in hypotheses}
        if missing:
            raise SystemExit(f"Unknown hypothesis key(s): {', '.join(sorted(missing))}")
    if not hypotheses:
        raise SystemExit("No hypotheses selected.")

    target_map = read_target_map(args.target_map)
    print(f"Loaded {len(hypotheses)} hypotheses; {len(target_map)} mapped targets.")
    print("Open Targets and STRING are public web APIs; use --dry-run to avoid API calls.")

    results = []
    for h in hypotheses:
        results.append(run_one(h, args, target_map))
    if not args.dry_run:
        write_index(results, args.out)

    print(f"\nAll outputs under {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
