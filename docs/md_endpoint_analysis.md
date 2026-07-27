# MD Endpoint Analysis

## Scope

The endpoint workflow operates on the raw corrected production directories. It
performs final-20-ns ligand-pose clustering for every corrected SGLT2 and OPRK1
trajectory. It prepares and runs energetic calculations only for the nine SGLT2
replicates:

- SGLT2/Mol_13144, replicates 1-3
- SGLT2/Mol_13733, replicates 1-3
- SGLT2/Mol_15088, replicates 1-3

The energetic comparison is restricted to these ligands because they share the
same receptor, membrane setup, protein force field and OpenFF parameterization
workflow. OPRK1/Mol_16614 is clustered, but its endpoint energy must not be
compared numerically with SGLT2.

## Environment

Clustering is a CPU analysis. AmberTools MM/PB(GB)SA also primarily uses CPU
cores; it is run on the GPU workstation because that machine holds the raw
trajectories and serialized systems.

Create a separate endpoint environment so AmberTools does not disturb the
CUDA-pinned production environment:

```bash
mamba env create -f environment-md-endpoint.yml
mamba activate anpdb-md-endpoint
MMPBSA.py --help >/dev/null
python -c "import mdtraj, openmm, parmed, scipy; print('endpoint imports OK')"
```

## Input Requirements

Every corrected run requires:

```text
md_runs/production/<system>/repN/
  config.toml
  equilibrated.pdb
  production.dcd
  production.log
  receptor_fixed.pdb
  system.xml
  system_solvated.pdb
```

`system.xml` is authoritative for covalent and nonbonded parameters. The script
fails rather than proceeding if ParmEd reports an unsupported OpenMM energy
term, if the atom counts differ, if the complex/receptor/ligand charges are
inconsistent, or if an omitted equilibration restraint has a nonzero default.

## Commands

Run every analysis currently possible and skip the two unfinished trajectories:

```bash
mamba activate anpdb-md-endpoint
python scripts/run_primary_md_endpoint_analysis.py \
  --scope completed \
  --mode run \
  --mpi-ranks 8
```

After all 12 corrected trajectories are complete:

```bash
python scripts/run_primary_md_endpoint_analysis.py \
  --scope all \
  --mode run \
  --mpi-ranks 8
```

Useful staged modes are:

```bash
# Clustering only; AmberTools is not invoked.
python scripts/run_primary_md_endpoint_analysis.py --scope completed --mode cluster

# Cluster and generate all Amber topologies, snapshots and MMPBSA.py inputs.
python scripts/run_primary_md_endpoint_analysis.py --scope completed --mode prepare

# Reuse cluster outputs and run/re-run the endpoint calculations.
python scripts/run_primary_md_endpoint_analysis.py \
  --scope completed \
  --mode run \
  --skip-clustering \
  --mpi-ranks 8
```

The per-run entry point is also available:

```bash
python scripts/md_endpoint_analysis.py \
  --run-dir md_runs/production/sglt2_mol13144/rep1 \
  --ligand-resname UNK \
  --prepare-mmpbsa \
  --run-mmpbsa \
  --mpi-ranks 8
```

## Clustering Method

The complete final 20 ns is streamed from the DCD in small chunks. Each selected
frame is periodically imaged around the receptor and aligned on the
transmembrane helical C-alpha core used by the primary MD analysis. Pairwise
ligand heavy-atom pose RMSDs are calculated without fitting the ligand to
itself. Average-linkage hierarchical clustering is cut at 2.0 A by default.

The dominant-cluster medoid is exported as the representative late-stage
structure. Protein residue numbering in medoid PDBs is restored from each
run's `receptor_fixed.pdb`; the mapping is recorded in
`residue_number_map.csv`.

## Endpoint-Energy Method

One hundred evenly spaced snapshots are selected from the final 20 ns of every
SGLT2 replicate. Snapshots are imaged, receptor-core aligned and translated so
that the explicit lipid-center z coordinate becomes zero before lipids, water
and ions are removed. ParmEd converts the exact OpenMM `system.xml`
parameterization into matched complex, receptor and ligand Amber topologies.
The `mbondi2` radii set is assigned and retained with `radiopt=0`.

MMPBSA.py performs two single-trajectory estimates:

- MM/GBSA with `igb=5` and 0.15 M salt as an aqueous sensitivity model.
- MM/PBSA with a 40 A implicit membrane slab centered at z=0, solute
  dielectric 4, membrane dielectric 7, solvent dielectric 80 and 0.15 M salt.

The implicit-membrane PB setup follows the Amber PBSA membrane model. The
underlying approach and single-trajectory approximation are described in the
[Amber MM/PBSA tutorial](https://ambermd.org/tutorials/advanced/tutorial3/index.php)
and the [MMPBSA.py paper](https://doi.org/10.1021/ct300418h).

These values are endpoint estimates without configurational entropy. They are
supportive comparative descriptors, not exact binding free energies,
experimental validation or evidence of pharmacological activity. The
independent replicate mean is the unit of replication.

## Outputs

Full working files remain under the ignored run directory:

```text
md_runs/production/<system>/repN/endpoint_analysis/
```

The batch script copies lightweight manuscript outputs to:

```text
results/md_analysis/manuscript_v2/endpoint_analysis/
```

This tracked directory contains cluster assignments and summaries, dominant
medoid PDBs, residue-number mappings, native MMPBSA.py results, per-replicate
endpoint estimates and across-replicate summaries. It deliberately excludes
raw DCDs, NetCDF snapshot trajectories and generated topology files.
