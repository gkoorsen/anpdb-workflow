# Amber-Prepared Production MD Configs

Use these configs for systems where cofactors are handled during curated system preparation instead of as separate OpenMM ffxml residue templates.

Current configs:

- `CYP1B1 / Mol_11315`: generated from AmberTools ff14SB/GAFF2 plus Shahrokh IC6 P450 heme parameters.
- `MAO-B / Mol_14056`: generated from AmberTools ff14SB/GAFF2 plus a CYF residue for oxidized FAD covalently linked to Cys397 through SG-C8M.

Generate the systems from the repo root with the AmberTools-capable `Docking` environment:

```bash
conda run -n Docking python scripts/md_build_amber_systems.py --target cyp1b1 --force
conda run -n Docking python scripts/md_build_amber_systems.py --target maob --force
```

The MAO-B builder creates `data/md_inputs/cofactors/amber/CYF_cys397_fad.mol2` and `.frcmod` from a capped Cys397-FAD AM1-BCC model. It does not use a free-FAD placeholder.

Expected Git LFS input bundle:

```text
data/md_inputs/amber_systems/
  cyp1b1_mol11315.prmtop
  cyp1b1_mol11315.inpcrd
  maob_mol14056.prmtop
  maob_mol14056.inpcrd
```

After generating those files, add and push them from the repo root:

```bash
git lfs install
git add data/md_inputs/amber_systems/cyp1b1_mol11315.prmtop data/md_inputs/amber_systems/cyp1b1_mol11315.inpcrd
git add data/md_inputs/amber_systems/maob_mol14056.prmtop data/md_inputs/amber_systems/maob_mol14056.inpcrd
git commit -m "Add curated Amber MD input systems"
git push
```

On the GPU machine:

```bash
git lfs install
git pull
git lfs pull
```

Dry-run first replicates:

```bash
python scripts/md_check_amber_inputs.py
python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --dry-run
python scripts/md_production_amber.py --config configs/md/amber_production/maob_mol14056_100ns_rep1.toml --dry-run
```

Run one replicate:

```bash
python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --equilibrate-only
python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --resume
python scripts/md_production_amber.py --config configs/md/amber_production/maob_mol14056_100ns_rep1.toml --equilibrate-only
python scripts/md_production_amber.py --config configs/md/amber_production/maob_mol14056_100ns_rep1.toml --resume
```

Equilibrate all ready Amber-prepared configs sequentially:

```bash
python scripts/md_batch_amber.py configs/md/amber_production/*.toml --equilibrate-only
```

After inspecting the equilibrated structures and manifests, run production sequentially:

```bash
python scripts/md_batch_amber.py configs/md/amber_production/*.toml --resume
```
