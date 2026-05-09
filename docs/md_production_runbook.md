# Production MD Runbook

## Local CPU Machine

Use this machine for source edits, syntax checks, and dry-runs only. Do not start 100 ns production jobs here.

```bash
python3 -m py_compile scripts/check_openmm_cuda.py scripts/md_prepare_inputs.py scripts/md_fetch_orient_sglt2.py scripts/md_check_inputs.py scripts/md_check_amber_inputs.py scripts/md_production.py scripts/md_production_amber.py scripts/md_batch.py scripts/md_batch_amber.py scripts/md_analyze_production.py
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
  amber_systems/
    cyp1b1_mol11315.prmtop
    cyp1b1_mol11315.inpcrd
    maob_mol14056.prmtop
    maob_mol14056.inpcrd
  cofactors/
    amber/
      CYF_cys397_fad.mol2
      CYF_cys397_fad.frcmod
      CYF_cys397_fad_manifest.json
```

The SGLT2 receptor and ligand pose must be in the same membrane-oriented coordinate frame before use. For CYP1B1 and MAO-B, prefer the Amber-prepared `amber_systems/` files. `cofactors/heme.xml` and `cofactors/fad.xml` are only needed if using the older OpenMM ffxml assembly configs under `configs/md/production/`.

Check readiness at any time:

```bash
python scripts/md_check_inputs.py
```

## Cofactor Systems

For CYP1B1 and MAO-B, the preferred route is full Amber-prepared systems instead of separate cofactor XML files. See `docs/cofactor_parameterization.md`.

Generate the curated CYP1B1 Amber system locally from the repo root:

```bash
conda run -n Docking python scripts/md_build_amber_systems.py --target cyp1b1 --force
```

The current CYP1B1 system uses Amber ff14SB, GAFF2/AM1-BCC ligand parameters, Shahrokh IC6 P450 heme parameters, TIP3P water, and approximately 0.15 M NaCl.

Generate the curated MAO-B Amber system locally from the repo root:

```bash
conda run -n Docking python scripts/md_build_amber_systems.py --target maob --force
```

The current MAO-B system uses Amber ff14SB, GAFF2/AM1-BCC ligand parameters, a tracked CYF Cys397-FAD residue generated from a capped AM1-BCC model, TIP3P water, and approximately 0.15 M NaCl. The repo intentionally does not create a free-FAD placeholder for production work.

Place curated Amber systems here:

```text
data/md_inputs/amber_systems/
  cyp1b1_mol11315.prmtop
  cyp1b1_mol11315.inpcrd
  maob_mol14056.prmtop
  maob_mol14056.inpcrd
```

These files are intentionally allowed into the repository through Git LFS. After generating or receiving the final curated systems, run:

```bash
git lfs install
git add data/md_inputs/amber_systems/cyp1b1_mol11315.prmtop data/md_inputs/amber_systems/cyp1b1_mol11315.inpcrd
git add data/md_inputs/amber_systems/maob_mol14056.prmtop data/md_inputs/amber_systems/maob_mol14056.inpcrd
git commit -m "Add curated Amber MD input systems"
git push
```

On the GPU machine, fetch the LFS payloads after cloning or pulling:

```bash
git lfs install
git lfs pull
```

Dry-run those configs:

```bash
python scripts/md_check_amber_inputs.py
python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --dry-run
python scripts/md_production_amber.py --config configs/md/amber_production/maob_mol14056_100ns_rep1.toml --dry-run
```

## Running

Dry-run the production plan:

```bash
python scripts/md_batch.py configs/md/production/*_rep1.toml --dry-run
```

Run one replicate:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --equilibrate-only
```

Inspect `equilibrated.pdb` and `equilibration_manifest.json`, then start production from the saved `production.chk`:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --resume
```

Equilibrate all first replicates sequentially on one GPU:

```bash
python scripts/md_batch.py configs/md/production/sglt2_*_rep1.toml --equilibrate-only
```

After inspection, run SGLT2 production sequentially:

```bash
python scripts/md_batch.py configs/md/production/sglt2_*.toml --resume
```

Equilibrate ready Amber-prepared CYP1B1 and MAO-B configs sequentially:

```bash
python scripts/md_batch_amber.py configs/md/amber_production/*.toml --equilibrate-only
```

After inspection, run Amber-prepared production sequentially:

```bash
python scripts/md_batch_amber.py configs/md/amber_production/*.toml --resume
```

## Outputs

Each run writes to `md_runs/production/<target_compound>/repN/`, including:

- `config.toml`
- `plan.json`
- `system_solvated.pdb`
- `minimized.pdb`
- `equilibrated.pdb`
- `equilibrated_state.xml`
- `equilibration_manifest.json`
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
