# Production MD Runbook

## Local CPU Machine

Use this machine for source edits, syntax checks, and dry-runs only. Do not start 100 ns production jobs here.

```bash
python3 -m py_compile scripts/check_openmm_cuda.py scripts/md_prepare_inputs.py scripts/md_fetch_orient_sglt2.py scripts/md_check_inputs.py scripts/md_production.py scripts/md_batch.py scripts/md_analyze_production.py
python3 scripts/md_prepare_inputs.py
python3 scripts/md_check_inputs.py
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

Stage ignored inputs under `data/md_inputs/` on the GPU machine. If you have copied the current local `output/` and `data/` folders onto that machine, start with:

```bash
python scripts/md_prepare_inputs.py
python scripts/md_fetch_orient_sglt2.py
python scripts/md_check_inputs.py
```

The preparer copies the local docking poses, SMILES table, CYP1B1 receptor with HEM retained, MAO-B receptor with FAD retained, and the unoriented 7VSI reference. The SGLT2 orientation script downloads the OPM/EncoMPASS membrane-oriented 7VSI coordinates, keeps chain A, removes MAP17 and OPM dummy atoms, and applies the same rigid-body transform to the docked SGLT2 ligand pose. Neither script fabricates cofactor force-field XML files.

The SGLT2 orientation script uses NumPy, so run it from the `anpdb-md` environment or another environment with NumPy installed.

Final expected layout:

```text
data/md_inputs/
  anpdb_truly_novel_std.csv
  poses/
    Mol_11315_4I8V_out.pdbqt
    Mol_13144_7VSI_out.pdbqt
    Mol_13144_7VSI_opm_oriented_out.pdbqt
    Mol_14056_2V5Z_out.pdbqt
  receptors/
    4I8V_chainA_heme_prepared.pdb
    7VSI_opm_oriented_clean.pdb
    2V5Z_chainA_fad_prepared.pdb
  cofactors/
    heme.xml
    fad.xml
```

The SGLT2 receptor and ligand pose must be in the same membrane-oriented coordinate frame before use. The CYP1B1 and MAO-B receptor files should retain their cofactors, and the cofactor XML files should match the residue names present in those PDB files. For the current MAO-B structure, the retained flavin cofactor is FAD.

Check readiness at any time:

```bash
python scripts/md_check_inputs.py
```

## Cofactor Parameter Files

The remaining CYP1B1 and MAO-B blockers are true parameterization inputs:

```text
data/md_inputs/cofactors/heme.xml
data/md_inputs/cofactors/fad.xml
```

These should be generated from curated Amber/CHARMM-compatible cofactor parameterization, not from placeholder XML. For CYP1B1, use a heme model appropriate for the CYP450 heme environment and make sure the OpenMM ffxml defines a `HEM` residue template. For MAO-B, use oxidized FAD with the intended charge/protonation state documented and make sure the OpenMM ffxml defines a `FAD` residue template. After adding either file, rerun `python scripts/md_check_inputs.py`.

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
