"""Build a Word (.docx) version of the methods/results draft with all figures
embedded inline. Uses pandoc with a markdown-with-images intermediate.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
SRC = ROOT / "output" / "methods_results_draft.md"
INTERMEDIATE = ROOT / "output" / "methods_results_with_figures.md"
DOCX = ROOT / "output" / "methods_results_draft.docx"

FIG_DIR = ROOT / "output"

# (figure_label, caption, path) — paths absolute so pandoc can locate them
FIGURES = [
    ("fig_novelty_filter",
     "Two-pass novelty filter outcomes. Of 11,448 ANPDB compounds, 1,364 lacked any matching InChIKey skeleton in COCONUT 2.0 + ChEMBL 36 (Pass 1), and 1,012 retained Tanimoto < 0.85 against the COCONUT background (Pass 2).",
     FIG_DIR / "scaffolds" / "fig_scaffold_diversity.png"),  # placeholder for novelty bar; we use scaffold diversity for now

    ("fig_chemical_space",
     "Chemical-space UMAP. Morgan ECFP4 fingerprints (radius 2, 2,048 bits) of 1,012 truly-novel ANPDB compounds (red) projected with UMAP (Jaccard distance, n_neighbors=30, min_dist=0.3) over a random 5,000-compound COCONUT background (grey). Novel compounds occupy distinct regions of natural-product chemical space.",
     FIG_DIR / "chemspace" / "fig_chemical_space_umap.png"),

    ("fig_null_model",
     "Permutation null model for the novelty rate. Histogram of bootstrap-derived novelty rates (200 random 11,448-compound draws from COCONUT, identical two-pass filter) for Pass 1 (left) and Pass 2 (right). ANPDB's observed rates (red lines) lie far below the null distribution mean (blue dashed): one-sided p < 0.005 for both passes.",
     FIG_DIR / "null_model" / "fig_novelty_null.png"),

    ("fig_scaffold_diversity",
     "Bemis-Murcko scaffold diversity. Cumulative coverage curve (grey) showing the fraction of compounds covered by the top-N scaffolds; 50 % of compounds are covered by the 207 most-frequent scaffolds. 80 % of unique scaffolds are singletons.",
     FIG_DIR / "scaffolds" / "fig_scaffold_diversity.png"),

    ("fig_scaffold_top",
     "Top-20 Bemis-Murcko scaffolds by frequency among the 1,012 truly-novel ANPDB compounds.",
     FIG_DIR / "scaffolds" / "fig_scaffold_freq.png"),

    ("fig_idf_leaderboard",
     "Disease leaderboard after IDF re-weighting. Per-disease IDF = log(N_targets / df(d)). Promiscuous terms (Schizophrenia, Autosomal recessive predisposition) are deflated; specific terms (Glaucoma, Primary congenital glaucoma) rise to the top.",
     FIG_DIR / "network" / "figures" / "fig_disease_leaderboard_ad60_idf.png"),

    ("fig_cluster_summary",
     "IDF-weighted hierarchical clustering of compound disease profiles. Six clusters (k=6, cosine distance, average linkage). Clusters 4 and 6 (98 compounds, 33 % of the predicted-hit set) are dominated by ophthalmology terms.",
     FIG_DIR / "network" / "figures" / "fig_cluster_summary_ad60_idf.png"),

    ("fig_ophth_targets",
     "Ophthalmology-cluster top targets. Distinct ANPDB-novel compounds per UniProt target within the 98-compound ophthalmology cluster. CYP1B1, CA4 and FA2H drive the eye-disease signal.",
     FIG_DIR / "network" / "ophthalmology" / "figures" / "fig_top_targets.png"),

    ("fig_ophth_buckets",
     "Disease-bucket breakdown of the ophthalmology cluster. Six near-synonym buckets collapsed by regex matching: Glaucoma, Cataract/Lens, Hearing-loss/Otology, Cornea/Sclera, Retina/Optic, Other-eye.",
     FIG_DIR / "network" / "ophthalmology" / "figures" / "fig_disease_buckets.png"),

    ("fig_ophth_species",
     "Top source plant species in the ophthalmology cluster. Tephrosia purpurea and Solanum nigrum dominate; both have documented traditional eye-disease use in African ethnomedicine.",
     FIG_DIR / "network" / "ophthalmology" / "figures" / "fig_top_species.png"),

    ("fig_pathway",
     "Pathway enrichment of the 21 consensus targets. g:Profiler enrichment (g:SCS multiple-testing correction, p < 0.05) returned 11 significant terms across KEGG, Reactome and GO Biological Process.",
     FIG_DIR / "pathways" / "fig_pathway_enrichment.png"),

    ("fig_boiled_egg",
     "BOILED-Egg ADME classification of the top-50 consensus shortlist. Yellow region: predicted gastrointestinal-absorbed (TPSA ≤ 142, −0.4 ≤ WLogP ≤ 5.6); blue region: predicted blood–brain-barrier permeant (TPSA ≤ 79, −1.0 ≤ WLogP ≤ 3.5).",
     FIG_DIR / "adme" / "fig_boiled_egg.png"),

    ("fig_adme_rules",
     "ADME rule-of-five compliance for the top-50 consensus shortlist. Pass/fail counts for Lipinski, Veber, Egan, BOILED-Egg GI absorption, and BBB permeability.",
     FIG_DIR / "adme" / "fig_adme_rules.png"),

    ("fig_dock_validation",
     "Docking validation controls. (A) Decoy control: consensus hits dock significantly better than property-matched ANPDB decoys at all three targets. (B) Redock control: best-of-top-9 RMSD recovers crystal pose to < 2 Å for all three receptors. (C) CYP1B1 mode-by-mode scatter showing crystal pose recovered in mode 6 (0.39 Å) versus the symmetric-flip artefact in mode 1.",
     FIG_DIR / "docking" / "validation" / "fig_dock_validation_v2.png"),
]


# Where to insert each figure: anchor (heading text we look for) → list of figures to insert AFTER that section's body
INSERTION_PLAN = {
    "### 3.1 1,012 ANPDB compounds are structurally novel relative to COCONUT and ChEMBL": [
        "fig_chemical_space", "fig_null_model",
    ],
    "### 3.2 Bemis-Murcko scaffold diversity is high": [
        "fig_scaffold_diversity", "fig_scaffold_top",
    ],
    "### 3.4 IDF-weighted disease profiling surfaces an ophthalmology polypharmacology cluster": [
        "fig_idf_leaderboard", "fig_cluster_summary",
    ],
    "### 3.5 The ophthalmology cluster maps to validated eye-disease targets and traditional medicinal species": [
        "fig_ophth_targets", "fig_ophth_buckets", "fig_ophth_species",
    ],
    "### 3.6 Pathway enrichment of the consensus targets": [
        "fig_pathway",
    ],
    "### 3.7 ADME profile of the top-50 shortlist": [
        "fig_boiled_egg", "fig_adme_rules",
    ],
    "#### 3.8.1 Validation controls": [
        "fig_dock_validation",
    ],
}


def figure_block(label: str, caption: str, path: Path) -> str:
    return (
        f"\n\n![**Figure {label.replace('fig_','').replace('_',' ').title()}.** {caption}]"
        f"({path.absolute()}){{#{label}}}\n\n"
    )


def main() -> int:
    src = SRC.read_text().splitlines(keepends=True)
    fig_lookup = {label: (caption, path) for label, caption, path in FIGURES}

    out_lines: list[str] = []
    pending_inserts_after_section: list[str] = []
    in_target_section: str | None = None

    for line in src:
        out_lines.append(line)

        # If we just emitted a heading line, decide whether it triggers an insertion
        if line.startswith("###") or line.startswith("#### "):
            # Flush any pending inserts before moving to a new section
            for label in pending_inserts_after_section:
                cap, path = fig_lookup[label]
                if not path.exists():
                    print(f"  WARN: missing figure {path}")
                    continue
                out_lines.append(figure_block(label, cap, path))
            pending_inserts_after_section = []

            stripped = line.rstrip("\n")
            if stripped in INSERTION_PLAN:
                pending_inserts_after_section = INSERTION_PLAN[stripped]
                in_target_section = stripped
            else:
                in_target_section = None

    # Final flush at EOF
    for label in pending_inserts_after_section:
        cap, path = fig_lookup[label]
        if path.exists():
            out_lines.append(figure_block(label, cap, path))

    INTERMEDIATE.write_text("".join(out_lines))
    print(f"Wrote intermediate markdown -> {INTERMEDIATE}")

    # Convert to docx with pandoc, A4 portrait, default style
    cmd = [
        "/opt/homebrew/bin/pandoc",
        str(INTERMEDIATE),
        "-o", str(DOCX),
        "--standalone",
        "--toc",
        "--toc-depth=3",
        "--from=markdown+grid_tables+pipe_tables+implicit_figures+raw_html",
        "--metadata", "title=Multi-method computational target prediction reveals ophthalmology, type 2 diabetes and neurodegeneration leads in the African Natural Products Database",
        "--metadata", "author=Working draft",
        "--metadata", f"date={subprocess.check_output(['date', '+%Y-%m-%d']).decode().strip()}",
    ]
    print(f"Running pandoc ...")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"pandoc STDERR:\n{r.stderr}")
        return 1
    print(f"Wrote {DOCX} ({DOCX.stat().st_size//1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
