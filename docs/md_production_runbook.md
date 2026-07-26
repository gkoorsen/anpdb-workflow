# Production MD Runbook

## Local CPU Machine

Use this machine for source edits, syntax checks, and dry-runs only. Do not start 100 ns production jobs here.

```bash
python3 -m py_compile scripts/check_openmm_cuda.py scripts/md_prepare_inputs.py scripts/md_fetch_orient_sglt2.py scripts/md_fetch_orient_oprk1.py scripts/md_check_inputs.py scripts/md_check_amber_inputs.py scripts/md_production.py scripts/md_production_amber.py scripts/md_batch.py scripts/md_batch_amber.py scripts/md_analyze_production.py
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

For the portable SGLT2 Mol_13144 replicate-2 RunPod workflow, including the
transfer-package builder and cross-GPU XML-state restart, see
[`docs/runpod_sglt2_rep2.md`](runpod_sglt2_rep2.md).
The sequential portable workflow for SGLT2 Mol_13733 replicates 1–3 is in
[`docs/runpod_sglt2_mol13733.md`](runpod_sglt2_mol13733.md).
The queued portable workflow for SGLT2 Mol_15088 replicates 2–3 is in
[`docs/runpod_sglt2_mol15088_reps2-3.md`](runpod_sglt2_mol15088_reps2-3.md).

## Input Bundle

Stage ignored inputs under `data/md_inputs/` on the GPU machine. If you have copied the current local `output/` and `data/` folders onto that machine, start with:

```bash
python scripts/md_prepare_inputs.py
python scripts/md_fetch_orient_sglt2.py
python scripts/md_fetch_orient_oprk1.py
python scripts/md_check_inputs.py
```

The preparer copies the local docking poses, SMILES table, CYP1B1 receptor with HEM retained, MAO-B receptor with FAD retained, the unoriented 7VSI reference, and the unoriented 4DJH OPRK1 reference. The SGLT2 orientation script downloads the OPM/EncoMPASS membrane-oriented 7VSI coordinates, keeps chain A, removes MAP17 and OPM dummy atoms, and applies the same rigid-body transform to the docked SGLT2 ligand pose. The OPRK1 orientation script downloads the EncoMPASS/OPM-oriented 4DJH coordinates and applies the same rigid-body transform to the Mol_16614 pose. Neither script fabricates cofactor force-field XML files.

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
    Mol_16614_4DJH_out.pdbqt
    Mol_16614_4DJH_opm_oriented_out.pdbqt
  receptors/
    4I8V_chainA_heme_prepared.pdb
    7VSI_opm_oriented_clean.pdb
    2V5Z_chainA_fad_prepared.pdb
    4DJH_OPRK1_clean_unoriented_reference.pdb
    4DJH_OPRK1_opm_oriented_clean.pdb
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

The SGLT2 and OPRK1 receptors and ligand poses must be in their respective membrane-oriented coordinate frames before use. Inspect `equilibrated.pdb` after the first OPRK1 equilibration because 4DJH contains an engineered T4 lysozyme fusion segment. For CYP1B1 and MAO-B, prefer the Amber-prepared `amber_systems/` files. `cofactors/heme.xml` and `cofactors/fad.xml` are only needed if using the older OpenMM ffxml assembly configs under `configs/md/production/`.

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

### Corrected primary membrane protocol

The 12 primary SGLT2/OPRK1 simulations use protocol
`membrane-endpoint-v2`: true restrained NVT with pressure coupling disabled,
restrained membrane NPT, unrestrained membrane NPT, and production with an
XY-isotropic/Z-free `MonteCarloMembraneBarostat` at zero surface tension.
The production DCD, solvated PDB, serialized OpenMM System, and
`endpoint_analysis_manifest.json` support a consistent single-trajectory
MM/GBSA or membrane-aware MM/PBSA workflow.

```bash
scripts/run_primary_membrane_productions.sh --preflight
scripts/run_primary_membrane_productions.sh --equilibrate --archive-obsolete
# Inspect the corrected equilibrated structures and manifests, then:
scripts/run_primary_membrane_productions.sh --produce
```

The launcher refuses old isotropic-barostat checkpoints. The archive option
moves them into `md_runs/obsolete_protocol/<timestamp>/` before rebuilding.
Pass `--python /path/to/env/bin/python` on non-WSL systems if needed.

For the corrected primary membrane reruns on RunPod/Linux or WSL, use the
pending-equilibration launcher. It runs SGLT2/Mol_13144, SGLT2/Mol_13733,
SGLT2/Mol_15088, and OPRK1/Mol_16614, three replicates each, skipping runs
that already have `equilibration_manifest.json` and `production.chk`:

```bash
scripts/run_pending_md_equilibrations.sh --dry-run-only
scripts/run_pending_md_equilibrations.sh --archive
```

If the MD environment is not the active Python, pass its Python explicitly:

```bash
scripts/run_pending_md_equilibrations.sh --python /path/to/env/bin/python --archive
```

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

Equilibrate the OPRK1 first replicate:

```bash
python scripts/md_production.py --config configs/md/production/oprk1_mol16614_100ns_rep1.toml --equilibrate-only
```

Inspect `md_runs/production/oprk1_mol16614/rep1/equilibrated.pdb`. If the lipid placement is acceptable, run all OPRK1 replicates:

```bash
python scripts/md_batch.py configs/md/production/oprk1_mol16614_100ns_rep*.toml --resume
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
- `system.xml`
- `integrator.xml`
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

Analyze only completed trajectories. OpenFF-prepared ligands are written as
`UNK`; the Amber-prepared MAO-B ligand is `LIG`.

```bash
python scripts/md_analyze_production.py \
  --run-dir md_runs/production/sglt2_mol13144/rep1 \
  --ligand-resname UNK
```

For the 12 corrected SGLT2 and OPRK1 runs on the GPU machine:

```bash
python scripts/run_primary_md_analysis.py \
  --scope corrected \
  --skip-publication
```

After all 15 primary trajectories, including MAO-B/Mol_14056, are available on
one machine:

```bash
python scripts/run_primary_md_analysis.py --scope all
```

The per-run analyzer performs PBC imaging before receptor alignment and writes:

- `alignment_core_atoms.csv`, documenting the target-specific alignment
  reference;
- `rmsd_timeseries.csv`, containing protein alignment-core RMSD, global backbone
  RMSD, ligand pose RMSD and ligand internal RMSD as separate quantities;
- `rmsf_ca.csv`;
- `contact_occupancy.csv`;
- `pose_retention_timeseries.csv`;
- `ligand_geometry_timeseries.csv` and `ligand_geometry_bonds.csv`;
- `ligand_protein_hbond_occupancy.csv`;
- `thermodynamic_timeseries.csv` when a production state log is available;
- `thermodynamic_qc_summary.json`;
- `analysis_summary.json`.

For lipid-containing SGLT2 and OPRK1 systems, the default `auto` alignment mode
uses frame-0 DSSP alpha-helical protein C-alpha atoms within +/-1.5 nm of the
lipid-centre z coordinate. MAO-B defaults to all protein backbone atoms. Use
`--alignment-selection` only for a predefined, documented override.

Ligand geometry and hydrogen-bond topology use parameterized bonds from
`system.xml`. Completed Amber runs without serialized systems fall back to the
configured `prmtop`. PDB bonds are retained only for a connectivity comparison,
not as the authoritative ligand graph.

The thermodynamic summary reports the configured 1 atm Monte Carlo barostat
target separately from observed state-log quantities. The saved production logs
do not provide an instantaneous-pressure time series.
