# Production MD Configs

These configs define three 100 ns replicates for each current top lead:

- `CYP1B1 / Mol_11315`
- `SGLT2 / Mol_13144`
- `MAO-B / Mol_14056`

They are intentionally guarded. Production configs expect an ignored local input bundle under `data/md_inputs/` and will fail before running if required receptors, ligand poses, SMILES, or cofactor parameter files are missing.

You can stage the easy local inputs from existing `output/` files with:

```bash
python scripts/md_prepare_inputs.py
python scripts/md_fetch_orient_sglt2.py
python scripts/md_check_inputs.py
```

Final expected ignored input bundle:

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

For the current MAO-B structure, the retained flavin cofactor is FAD. If a future MAO target uses FMN instead, create a separate receptor/config pair and update the expected residue and parameter file accordingly.

The SGLT2 orientation script creates both `7VSI_opm_oriented_clean.pdb` and the matching transformed ligand pose. The heme and FAD XML files must come from curated cofactor parameterization and are intentionally not generated as placeholders.

Run one config:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml
```

Resume after interruption:

```bash
python scripts/md_production.py --config configs/md/production/sglt2_mol13144_100ns_rep1.toml --resume
```

Run a batch sequentially on one GPU:

```bash
python scripts/md_batch.py configs/md/production/*_rep1.toml --resume
```
