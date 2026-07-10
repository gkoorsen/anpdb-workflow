# MD Ligand Preparation Diagnosis

## Summary

The ligand-connectivity problem is not a LigPlot-only issue and is not a simple final-PDB export artifact. The root cause is an atom-mapping bug in the OpenMM/OpenFF ligand-preparation path used by `scripts/md_production.py`.

`scripts/md_production.py` calls `pose_to_openff_mol()` from `scripts/md_run.py`. In the previous implementation, RDKit's `GetSubstructMatch()` result was interpreted in the wrong direction when docked pose coordinates were grafted onto the canonical ligand molecule. This scrambled heavy-atom coordinates before OpenFF parameterisation and system construction.

The script has been patched so that:

- canonical-to-pose atom mapping is generated in the correct direction;
- the fallback reverse mapping is explicitly inverted;
- ligand heavy-heavy bond lengths are checked after coordinate grafting;
- future runs fail immediately if grafted ligand geometry is invalid.

## Evidence

`results/md_analysis/ligand_graft_mapping_diagnostics.csv` compares the old and fixed mapping paths on the MD input PDBQT poses.

Using the old mapping, the grafted ligands have many impossible heavy-atom bond lengths:

| Compound | Bad bonds with old mapping | Maximum bond length |
|---|---:|---:|
| Mol_13144 | 32 | 11.68 A |
| Mol_13733 | 29 | 8.82 A |
| Mol_15088 | 32 | 11.17 A |
| Mol_16614 | 24 | 9.44 A |

Using the fixed mapping, all tested OpenFF input poses have zero bad bonds.

`results/md_analysis/ligplot_fig6_fixed/ligand_geometry_diagnostics.csv` checks the existing MD endpoint PDBs. The endpoint ligand geometries are invalid for:

| Complex | Bad endpoint bonds | Maximum endpoint bond length |
|---|---:|---:|
| SGLT2 / Mol_13144 | 7 | 3.00 A |
| SGLT2 / Mol_13733 | 5 | 2.72 A |

The endpoint geometries for SGLT2 / Mol_15088 and OPRK1 / Mol_16614 are chemically contiguous, but those systems still passed through the affected OpenFF coordinate-grafting code and should be treated cautiously if the claim is specifically that MD validated the original docked pose.

MAO-B / Mol_14056 used the AmberTools `prmtop`/`inpcrd` workflow (`md_production_amber.py`) and is not affected by this OpenFF grafting bug.

## Recommendation

Do not use the existing SGLT2 / Mol_13144 and SGLT2 / Mol_13733 MD simulations for manuscript pose-stability claims.

For a clean manuscript MD set, rerun the OpenFF-prepared systems from the patched code:

- SGLT2 / Mol_13144
- SGLT2 / Mol_13733
- SGLT2 / Mol_15088
- OPRK1 / Mol_16614

MAO-B / Mol_14056 does not need rerunning for this reason.

