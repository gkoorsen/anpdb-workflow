# RunPod: SGLT2 Mol_15088 Replicates 2–3

This package queues the two already-equilibrated Mol_15088 replicas that are
not running locally. Replicate 2 runs first; replicate 3 starts automatically
after rep2 completes. The launcher skips completed runs and resumes an
interrupted run from its RunPod-generated checkpoint.

## Install the package on the existing pod

Do this after the preceding RunPod queue has completed:

```bash
cd /workspace/anpdb-workflow
git pull --ff-only origin agent/runpod-sglt2-rep2

sha256sum -c \
  transfer/runpod/sglt2_mol15088_reps2-3_portable_*.tar.gz.sha256
tar -xzf transfer/runpod/sglt2_mol15088_reps2-3_portable_*.tar.gz
```

## Preflight both replicas

```bash
PYTHON=/workspace/miniforge3/envs/anpdb-md/bin/python \
  scripts/run_runpod_sglt2_mol15088_reps2-3.sh --check-only
```

The output must report `CUDA context: OK`, `input_ready: true` for both
configs, and the final preflight-passed message.

## Add the queue in tmux

Do not run this concurrently with another MD on the same GPU.

```bash
tmux new-session -s sglt2_mol15088
cd /workspace/anpdb-workflow
PYTHON=/workspace/miniforge3/envs/anpdb-md/bin/python \
  scripts/run_runpod_sglt2_mol15088_reps2-3.sh
```

Detach with `Ctrl-b`, then `d`. Replicate 3 starts automatically after
replicate 2 completes.

If only SSH disconnects, reattach with:

```bash
tmux attach -t sglt2_mol15088
```

If the pod process stops, rerun the same launcher in a new tmux session. It
skips completed work and resumes partial work only when a RunPod checkpoint is
present. Do not alter a partial production directory.

## Retrieve final outputs

After both logs end at step `53750000`:

```bash
tar -czf /workspace/sglt2_mol15088_reps2-3_completed.tar.gz \
  md_runs/production/sglt2_mol15088/rep2 \
  md_runs/production/sglt2_mol15088/rep3
sha256sum /workspace/sglt2_mol15088_reps2-3_completed.tar.gz \
  > /workspace/sglt2_mol15088_reps2-3_completed.tar.gz.sha256
```

Download and checksum-verify the archive before terminating the pod.
