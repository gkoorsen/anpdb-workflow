# Disease-Network Pharmacology Branch

This branch is a disease-specific follow-up to the existing ANPDB novelty and consensus target-prediction workflow. It is intentionally parallel to the MD work and does not change the docking or MD inputs.

## Aim

For each lead hypothesis:

- CYP1B1 / glaucoma
- SGLT2 / type 2 diabetes mellitus
- MAO-B / Parkinson disease

the branch asks whether the predicted drug target sits inside, overlaps, or is close to an external disease-target network.

## Data Sources

- Open Targets Platform GraphQL API: disease-associated targets by Open Targets disease identifier.
- STRING API: human protein-protein interaction network for the disease module plus predicted drug target.
- Existing local consensus set: `output/chembl_nn/consensus_pidgin_chembl_nn_ad60.tsv`.

The disease IDs and parameters are configured in:

- `configs/network_pharmacology/lead_hypotheses.csv`
- `configs/network_pharmacology/target_map.csv`

## Run

```bash
python scripts/network_pharmacology_branch.py
```

Useful options:

```bash
python scripts/network_pharmacology_branch.py --force
python scripts/network_pharmacology_branch.py --skip-enrichment
python scripts/network_pharmacology_branch.py --hypothesis sglt2_t2d
```

## Outputs

Each hypothesis writes to `output/network_pharmacology/<hypothesis>/`:

- `disease_targets_opentargets.csv`: retained disease targets from Open Targets.
- `drug_targets.csv`: configured lead target plus any mapped consensus co-targets.
- `ppi_edges_string.tsv`: STRING PPI edges among disease targets and drug targets.
- `ppi_nodes_ranked.csv`: node-level disease/drug annotations and graph metrics.
- `hub_gene_candidates.csv`: top-ranked hub/context genes.
- `network_proximity.csv`: observed drug-target-to-disease-module proximity and empirical randomization p-value.
- `enrichment_results.csv`: optional g:Profiler enrichment of overlap and nearest-neighbour module genes.
- `summary.md`: short manuscript-facing interpretation.

## Interpretation Rules

Treat this branch as mechanistic context, not proof of activity.

- A direct overlap means the predicted target is already an Open Targets disease-associated target.
- A short PPI distance means the predicted target is near the disease module in STRING.
- A hub gene is a network-context gene, not necessarily a compound target.
- Claims still need the MD/biochemical validation branch before being phrased as mechanism of action.
