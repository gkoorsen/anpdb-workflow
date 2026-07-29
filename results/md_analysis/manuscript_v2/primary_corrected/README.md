# Corrected primary MD analysis

This directory contains the validated, lightweight per-replicate outputs for
the completed corrected membrane trajectories:

- SGLT2/Mol_13144, replicates 1–3
- SGLT2/Mol_13733, replicates 1–3
- SGLT2/Mol_15088, replicates 1–3
- OPRK1/Mol_16614, replicates 1–3

Each trajectory contains 2,000 frames spanning approximately 100 ns. Analysis
used all frames, parameterized connectivity from `system.xml`, periodic imaging
around the receptor, and membrane-core C-alpha alignment. The analyzer streamed
25 frames at a time to remain safe on a 16 GB workstation.

Each replicate includes full-resolution RMSD, pose-retention, ligand-geometry
and thermodynamic tables; C-alpha RMSF; overall and 20 ns block contact
occupancies; hydrogen-bond occupancies; five exact 400-frame block summaries;
and a 400-frame final-20-ns summary. Raw trajectories and serialized systems are
deliberately excluded from Git.

The corrected primary set now contains all 12 completed membrane trajectories.
