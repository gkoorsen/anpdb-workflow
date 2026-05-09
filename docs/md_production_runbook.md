# Production MD Runbook

## Local CPU Machine

Use this machine for source edits, syntax checks, and dry-runs only. Do not start 100 ns production jobs here.

```bash
python3 -m py_compile scripts/check_openmm_cuda.py scripts/md_production.py scripts/md_batch.py scripts/md_analyze_production.py
python3 scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --dry-run
python3 scripts/md_batch.py configs/md/production/*_rep1.toml --dry-run
```

The dev config can be used for a tiny apo/water smoke test if the full MD environment is installed:

```bash
python scripts/md_production.py --config configs/md/dev/apo_water_10ps.toml --allow-non-production
```

## GPU WSL Setup

Clone the repo:

```bash
git clone https://github.com/gkoorsen/anpdb-workflow.git
cd anpdb-workflow
```

Create the environment:

```bash
mamba env create -f environment-md.yml
mamba activate anpdb-md
python scripts/check_openmm_cuda.py --device-index 0 --precision mixed
```

If CUDA is unavailable, fix the NVIDIA driver/WSL/CUDA/OpenMM stack before running production configs.

## Input Bundle

Stage ignored inputs under `data/md_inputs/` on the GPU machine:

```text
data/md_inputs/
  anpdb_truly_novel_std.csv
  poses/
    Mol_11315_4I8V_out.pdbqt
    Mol_13144_7VSI_out.pdbqt
    Mol_14056_2V5Z_out.pdbqt
  receptors/
    4I8V_chainA_heme_prepared.pdb
    7VSI_opm_oriented_clean.pdb
    2V5Z_chainA_fad_prepared.pdb
  cofactors/
    heme.xml
    fad.xml
```

The SGLT2 receptor should be membrane-oriented before use. The CYP1B1 and MAO-B receptor files should retain their cofactors, and the cofactor XML files should match the residue names present in those PDB files.

## Running

Dry-run the production plan:

```bash
python scripts/md_batch.py configs/md/production/*_rep1.toml --dry-run
```

Run one replicate:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml
```

Resume safely after interruption:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --resume
```

Run all first replicates sequentially on one GPU:

```bash
python scripts/md_batch.py configs/md/production/*_rep1.toml --resume
```

Run all configs sequentially:

```bash
python scripts/md_batch.py configs/md/production/*.toml --resume
```

## Outputs

Each run writes to `md_runs/production/<target_compound>/repN/`, including:

- `config.toml`
- `plan.json`
- `system_solvated.pdb`
- `minimized.pdb`
- `equilibrated.pdb`
- `production.dcd`
- `production.log`
- `production.chk`
- `production_plan.json`
- `final.pdb`
- `final_state.xml`
- `run_manifest.json`

These files are ignored by Git.

## Analysis

After a run completes:

```bash
python scripts/md_analyze_production.py --run-dir md_runs/production/sglt2_mol13144/rep1 --ligand-resname UNL
```

This writes RMSD, C-alpha RMSF, contact occupancy, and a JSON summary under the run's `analysis/` folder.
