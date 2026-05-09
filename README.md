# ANPDB Novelty

Python workflows for ANPDB novelty filtering, target-prioritisation, docking, and molecular-dynamics follow-up.

## Repository Scope

This repository is intended to track source code, lightweight configuration, runbooks, and curated publication-ready summaries. Large external datasets, generated outputs, docking intermediates, and MD trajectories are intentionally excluded from Git.

Ignored by default:

- `data/`: raw ANPDB/ChEMBL/COCONUT inputs and derived fingerprint/index files.
- `output/`: generated tables, figures, docking files, and previous MD artifacts.
- `runs/` and `md_runs/`: long production MD outputs, checkpoints, and trajectories.
- `ppt/node_modules/`: local JavaScript dependencies.

For publication and reproducibility, large artifacts should be archived separately with checksums, for example through Zenodo, OSF, institutional storage, or a project data bucket.

## Current MD Direction

The short OpenMM MD workflow in `scripts/md_run.py` was built for initial pose-stability checks. The next production workflow should use explicit 100 ns GPU runs with checkpointing, run manifests, and target-specific biological setup:

- `SGLT2`: membrane-embedded simulation.
- `CYP1B1`: publication-grade runs should retain and parameterise heme.
- `MAO-B`: publication-grade runs should retain and parameterise FAD.

See `docs/repo_bootstrap.md` for GitHub and WSL setup notes.
