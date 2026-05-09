# Amber-Prepared Production MD Configs

Use these configs for systems where cofactors are handled during curated system preparation instead of as separate OpenMM ffxml residue templates.

Current configs:

- `CYP1B1 / Mol_11315`: expects heme and Cys457-Fe coordination in the prepared Amber system.
- `MAO-B / Mol_14056`: expects oxidized FAD covalently patched to Cys397 in the prepared Amber system.

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
python scripts/md_batch_amber.py configs/md/amber_production/*_rep1.toml --dry-run
```

Run one replicate:

```bash
python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml
```

Run all Amber-prepared configs sequentially:

```bash
python scripts/md_batch_amber.py configs/md/amber_production/*.toml --resume
```
