# RunPod: SGLT2 Mol_13144 Replicate 2

This procedure starts the already-equilibrated `membrane-endpoint-v2` system
on RunPod and produces the same 100 ns outputs needed for trajectory and
MM/GBSA or membrane-aware MM/PBSA analysis.

The portable package deliberately uses `equilibrated_state.xml`, not the local
binary `production.chk`. OpenMM XML state includes the step count, positions,
velocities, periodic box and context parameters and is portable across CUDA
GPU models. A binary checkpoint is best reserved for resuming on the same
software/hardware stack after a later interruption.

## 1. Create the transfer package locally

From the repository root:

```bash
chmod +x scripts/package_runpod_sglt2_rep2.sh scripts/run_runpod_sglt2_rep2.sh
scripts/package_runpod_sglt2_rep2.sh
```

This creates an archive and matching SHA-256 file under `transfer/runpod/`.
The archive contains the exact serialized OpenMM System and integrator,
portable equilibrated state, topology, analysis selections, config, required
input metadata, environment definition, runner and this guide. `md_runs/` is
Git-ignored, so this archive must be uploaded separately; a Git clone alone is
not sufficient.

## 2. Create the RunPod pod

Use a CUDA-capable NVIDIA GPU pod with a persistent volume mounted at
`/workspace`. The workload needs roughly 8 GB GPU memory or more and at least
25 GB persistent disk for the environment, package, checkpoints and the final
trajectory. More disk is preferable. Start a terminal in the pod and confirm:

```bash
nvidia-smi
df -h /workspace
```

Upload both the `.tar.gz` and `.tar.gz.sha256` files into `/workspace` using
the RunPod file uploader, SCP/SSH, or your normal object-storage transfer.

## 3. Verify and extract

Replace the archive name below with the uploaded filename:

```bash
cd /workspace
sha256sum -c sglt2_mol13144_rep2_portable_*.tar.gz.sha256
mkdir -p anpdb-workflow
tar -xzf sglt2_mol13144_rep2_portable_*.tar.gz -C anpdb-workflow
cd anpdb-workflow
```

Do not extract this package over a directory containing an existing rep2
`production.dcd` or `production.log`.

## 4. Install the exact MD environment

If the pod does not already provide micromamba, install Miniforge or
micromamba, then create the environment:

```bash
micromamba create -y -n anpdb-md -f environment-md.yml
micromamba activate anpdb-md
```

If shell activation is unavailable, prefix commands with:

```bash
micromamba run -n anpdb-md <command>
```

Verify imports, CUDA, all required files and the config without starting MD:

```bash
scripts/run_runpod_sglt2_rep2.sh --check-only
```

The preflight must print `CUDA context: OK`. Do not proceed if it uses the CPU
or Reference platform.

## 5. Start production in tmux

```bash
tmux new-session -s sglt2_rep2
scripts/run_runpod_sglt2_rep2.sh
```

Detach with `Ctrl-b`, then `d`. Reattach with:

```bash
tmux attach -t sglt2_rep2
```

Monitor progress in a second terminal:

```bash
tail -f md_runs/production/sglt2_mol13144/rep2/production.log
nvidia-smi
```

The first production record appears after 50 ps. On completion the run
directory must contain `production.dcd`, `production.log`, `production.chk`,
`final.pdb`, `final_state.xml` and `run_manifest.json`.

## 6. Resume after an interruption on the same pod

Once production has started, resume from the new RunPod-generated binary
checkpoint—not from the equilibration XML:

```bash
python scripts/md_production.py \
  --config configs/md/production/sglt2_mol13144_100ns_rep2.toml \
  --resume
```

Never use `--resume-state` after `production.dcd` or `production.log` exists;
the runner refuses this to prevent overwriting or duplicating frames.

## 7. Retrieve and verify outputs

Copy the complete directory back before terminating the pod:

```text
md_runs/production/sglt2_mol13144/rep2/
```

At minimum retain the trajectory, log, latest checkpoint, final PDB/state,
run manifest, serialized System/integrator, solvated topology, production plan
and endpoint-analysis manifest. Generate a checksum before transfer:

```bash
tar -czf /workspace/sglt2_mol13144_rep2_completed.tar.gz \
  md_runs/production/sglt2_mol13144/rep2
sha256sum /workspace/sglt2_mol13144_rep2_completed.tar.gz \
  > /workspace/sglt2_mol13144_rep2_completed.tar.gz.sha256
```

Verify the checksum after downloading. Do not terminate the pod until the
downloaded archive has been tested successfully.
