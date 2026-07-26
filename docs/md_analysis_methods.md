# Molecular Dynamics Analysis Methods

## Primary Manuscript Set

The primary manuscript analysis comprises three independent 100 ns production
replicates for five complexes: SGLT2/Mol_13144, SGLT2/Mol_13733,
SGLT2/Mol_15088, MAO-B/Mol_14056 and OPRK1/Mol_16614. The 12 SGLT2 and OPRK1
trajectories must come from the corrected OpenFF ligand-coordinate grafting
workflow. MAO-B/Mol_14056 used the unaffected Amber-prepared system. CYP1B1 is
excluded from this comparative analysis.

## Run and Thermodynamic Quality Control

Run completeness is checked from the production manifest, trajectory, state log
and final structure. The production log is summarized for every available
thermodynamic field, including temperature, density, box volume and potential,
kinetic and total energy when those fields were recorded. Means, standard
deviations, ranges and linear slopes over time are retained as quality-control
descriptors.

The simulations used a Monte Carlo barostat configured with a target pressure of
1 atm. The production logs do not contain instantaneous pressure, and a pressure
time series cannot be reconstructed from the saved DCD coordinates. Accordingly,
1 atm is reported only as the configured barostat target. Stability assessment
uses the observed density and volume traces and does not claim that instantaneous
pressure was measured.

Ligand covalent geometry is checked in the equilibrated, intermediate and final
coordinates against the parameterized topology before structural analysis.
Membrane systems are inspected to confirm intact protein-ligand complexes,
appropriate membrane placement and retention of the ligand in the receptor
binding region.

## Trajectory Preprocessing

Trajectories are analysed with MDTraj, pandas and NumPy. Before alignment or
distance calculations, periodic trajectories are imaged around the protein so
that the receptor and ligand occupy the same periodic image. Each imaged
trajectory is then superposed on production frame 0 using protein backbone atoms.
The original DCD is never modified.

All frames are used for time-series plots. At the configured 50 ps trajectory
interval, a complete 100 ns trajectory contains approximately 2,000 frames.
Final-20-ns summaries are used for late-stage pose and contact retention, and
20 ns block summaries are retained to assess convergence without treating
individual frames as independent replicates.

## Structural Stability and Pose Retention

Protein backbone RMSD is calculated directly from the receptor-aligned
coordinates. Per-residue C-alpha RMSF is calculated around each atom's mean
position after receptor alignment.

Two different ligand RMSD quantities are retained and must not be conflated:

- **Ligand pose RMSD** is the direct heavy-atom RMSD from production frame 0
  after receptor alignment, without fitting the ligand to itself. This measures
  movement and conformational change within the binding site.
- **Ligand internal RMSD** is calculated after fitting ligand heavy atoms to the
  ligand reference. This measures internal conformational change only and is not
  described as binding-pose stability.

Additional pose-retention descriptors include ligand centre-of-mass displacement,
ligand-to-initial-pocket centre-of-mass distance, minimum ligand-protein
heavy-atom distance, ligand position relative to the membrane centre, contact
counts at 4, 6 and 8 A, and the fraction of the initial 4 A contacting residues
retained in each frame.

## Contacts and Hydrogen Bonds

Residue contact occupancy is the fraction of analysed frames in which at least
one residue heavy atom lies within 4 A of a ligand heavy atom. When occupancies
are averaged across replicates, a residue absent from a replicate contact table
contributes zero rather than being omitted. A consensus persistent contact can
be defined as occupancy of at least 50% in at least two of three replicates.

Candidate ligand-protein hydrogen-bond atom triplets are identified from the
bonded topology and evaluated frame by frame. A hydrogen bond is present when the
donor-acceptor distance is no greater than 3.5 A and the donor-hydrogen-acceptor
angle is at least 120 degrees. Occupancy is reported for each donor-acceptor
pair; hydrogen bonds with occupancy of at least 20% are highlighted.

## Replicate Summaries and Figures

The independent trajectory is the unit of replication. Complex-level values are
reported as the mean and standard deviation of the three replicate means.
Individual MD frames are not treated as independent observations for hypothesis
testing.

The main MD stability figure reports protein backbone RMSD, ligand pose RMSD,
late-stage pose/contact retention and consensus persistent contacts. Per-replicate
traces, C-alpha RMSF, ligand SASA, membrane-position quality control,
thermodynamic summaries, hydrogen-bond occupancies and 20 ns block summaries are
retained for supplementary reporting.

Representative interaction structures should be selected as medoids of dominant
late-stage ligand-pose clusters rather than as arbitrary final frames. Interaction
diagrams must use the parameterized ligand topology for connectivity and may be
annotated with trajectory-derived contact and hydrogen-bond occupancies.

## Energetic Interpretation

Direct ligand-protein Coulomb and Lennard-Jones interaction energies may be used
as within-system trajectory descriptors, but they are not binding free energies
because solvent, membrane, reorganization and entropic terms are omitted.
Absolute energetic values must not be compared between the OpenMM/OpenFF
membrane systems and the Amber-prepared MAO-B system.

If endpoint MM/GBSA calculations are added, the principal ranking comparison
should be restricted to the three SGLT2 ligands, which share a receptor and
parameterization workflow. Each replicate should be analysed independently and
replicate means and standard deviations reported. MM/GBSA estimates are
supportive computational evidence, not experimental validation.
