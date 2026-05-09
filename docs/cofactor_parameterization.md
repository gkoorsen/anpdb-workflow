# Cofactor Parameterization Notes

This repository now supports two production MD paths:

- OpenMM ffxml assembly for systems that can be built safely from receptor, ligand, and force-field files.
- Amber-prepared full systems for cofactors that need target-specific chemistry.

Use the Amber-prepared path for CYP1B1 and MAO-B unless you already have validated OpenMM ffxml templates for the exact cofactor states.

## CYP1B1 Heme

The 4I8V structure contains heme coordinated by Cys457 SG to HEM FE. Do not parameterize this as a generic free heme ligand.

Recommended preparation:

1. Start from the staged receptor and docked ligand:

   ```text
   data/md_inputs/receptors/4I8V_chainA_heme_prepared.pdb
   data/md_inputs/poses/Mol_11315_4I8V_out.pdbqt
   ```

2. Build the protein-ligand model in AmberTools or a trusted molecular editor, keeping HEM and Cys457 in the intended coordination state.

3. Use AmberTools MCPB.py for the metal center. The model should include the heme iron, porphyrin nitrogens, and the Cys457 sulfur coordination. Document the oxidation/spin/ligation state in the run notes.

4. Solvate and ionize the complete protein-ligand-cofactor system with `tleap`.

5. Save the final files as:

   ```text
   data/md_inputs/amber_systems/cyp1b1_mol11315.prmtop
   data/md_inputs/amber_systems/cyp1b1_mol11315.inpcrd
   ```

6. Add them with Git LFS:

   ```bash
   git lfs install
   git add data/md_inputs/amber_systems/cyp1b1_mol11315.prmtop data/md_inputs/amber_systems/cyp1b1_mol11315.inpcrd
   ```

7. Dry-run the Amber runner:

   ```bash
   python scripts/md_check_amber_inputs.py
   python scripts/md_production_amber.py --config configs/md/amber_production/cyp1b1_mol11315_100ns_rep1.toml --dry-run
   ```

## MAO-B FAD

The 2V5Z MAO-B structure contains FAD, and human MAO-B FAD is covalently attached to Cys397. Do not parameterize FAD as a free noncovalent ligand unless the scientific question explicitly calls for an apo/noncovalent control.

Recommended preparation:

1. Start from:

   ```text
   data/md_inputs/receptors/2V5Z_chainA_fad_prepared.pdb
   data/md_inputs/poses/Mol_14056_2V5Z_out.pdbqt
   ```

2. Build a covalent Cys397-FAD model using a curated Amber/CHARMM-compatible patch or residue definition. Document the FAD redox/protonation state.

3. Solvate and ionize the complete protein-ligand-cofactor system.

4. Save the final files as:

   ```text
   data/md_inputs/amber_systems/maob_mol14056.prmtop
   data/md_inputs/amber_systems/maob_mol14056.inpcrd
   ```

5. Add them with Git LFS:

   ```bash
   git lfs install
   git add data/md_inputs/amber_systems/maob_mol14056.prmtop data/md_inputs/amber_systems/maob_mol14056.inpcrd
   ```

6. Dry-run the Amber runner:

   ```bash
   python scripts/md_check_amber_inputs.py
   python scripts/md_production_amber.py --config configs/md/amber_production/maob_mol14056_100ns_rep1.toml --dry-run
   ```

## Why Full Amber Systems

OpenMM can read Amber `prmtop`/`inpcrd` directly. For heme iron coordination and covalent FAD, full-system Amber preparation avoids forcing complex cofactor chemistry into separate ffxml files after the fact. The production runner records the exact topology and coordinate paths in `run_manifest.json`.
