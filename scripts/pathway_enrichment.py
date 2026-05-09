"""KEGG / Reactome / GO pathway enrichment on consensus target set.

Uses g:Profiler REST API (no authentication needed) to run enrichment
on the 21 consensus UniProt targets from PIDGIN + ChEMBL-NN.

Outputs (under output/pathways/)
--------------------------------
- enrichment_results.csv       full g:Profiler results (KEGG, REAC, GO:BP)
- enrichment_summary.txt       top pathways per source
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from urllib import request, error

import pandas as pd

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
CONSENSUS = ROOT / "output" / "chembl_nn" / "consensus_pidgin_chembl_nn_ad60.tsv"
OUT = ROOT / "output" / "pathways"
OUT.mkdir(parents=True, exist_ok=True)

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def run_gprofiler(gene_list: list[str]) -> list[dict]:
    payload = {
        "organism": "hsapiens",
        "query": gene_list,
        "sources": ["KEGG", "REAC", "GO:BP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }
    data = json.dumps(payload).encode("utf-8")
    req = request.Request(
        GPROFILER_URL,
        data=data,
        headers={"Content-Type": "application/json"},
    )
    try:
        with request.urlopen(req, timeout=60) as resp:
            result = json.loads(resp.read().decode("utf-8"))
    except error.URLError as e:
        print(f"g:Profiler API error: {e}", file=sys.stderr)
        return []
    return result.get("result", [])


def main() -> int:
    consensus = pd.read_csv(CONSENSUS, sep="\t", dtype=str, keep_default_na=False)
    targets = sorted(consensus["uniprot"].unique())
    print(f"Consensus targets: {len(targets)}")
    print(f"  {', '.join(targets[:10])}{'...' if len(targets) > 10 else ''}")

    print(f"\nQuerying g:Profiler (KEGG + Reactome + GO:BP) ...")
    results = run_gprofiler(targets)
    print(f"  returned {len(results)} enriched terms")

    if not results:
        print("No enrichment results — writing empty output", file=sys.stderr)
        with open(OUT / "enrichment_summary.txt", "w") as fh:
            fh.write("No significant enrichment found.\n")
        return 0

    rows = []
    for r in results:
        rows.append({
            "source": r.get("source", ""),
            "term_id": r.get("native", ""),
            "term_name": r.get("name", ""),
            "p_value": r.get("p_value", 1.0),
            "term_size": r.get("term_size", 0),
            "intersection_size": r.get("intersection_size", 0),
            "query_size": r.get("query_size", 0),
            "intersections": " | ".join(
                str(x) if isinstance(x, str) else ", ".join(x) if isinstance(x, list) else str(x)
                for x in r.get("intersections", [])
            ),
        })

    df = pd.DataFrame(rows).sort_values(["source", "p_value"])
    df.to_csv(OUT / "enrichment_results.csv", index=False)

    summary = "Pathway Enrichment — Consensus Targets\n"
    summary += "=" * 40 + "\n\n"
    summary += f"Query: {len(targets)} UniProt IDs from PIDGIN + ChEMBL-NN consensus\n"
    summary += f"Method: g:Profiler (g:SCS multiple testing correction, p < 0.05)\n\n"

    for source in ["KEGG", "REAC", "GO:BP"]:
        sub = df[df["source"] == source]
        summary += f"\n{source} — {len(sub)} enriched terms\n"
        summary += "-" * 40 + "\n"
        for _, r in sub.head(15).iterrows():
            summary += (
                f"  {r['term_id']:<20} {r['term_name'][:55]:<55}  "
                f"p={r['p_value']:.2e}  "
                f"hits={r['intersection_size']}/{r['term_size']}\n"
            )

    with open(OUT / "enrichment_summary.txt", "w") as fh:
        fh.write(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
