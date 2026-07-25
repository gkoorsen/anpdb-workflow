# RunPod: SGLT2 Mol_13733 Replicates 1–3

This package runs all three already-equilibrated `membrane-endpoint-v2`
replicates sequentially on one RunPod GPU. Untouched replicates start from
portable XML states. Interrupted replicates resume from checkpoints generated
on RunPod, and completed replicates are skipped.

## After the Mol_13144 replicate-2 run finishes

First verify that its last log step is `53750000` and that `final_state.xml`
and `run_manifest.json` exist. Archive it before moving on:

```bash
cd /workspace/anpdb-workflow
tar -czf /workspace/sglt2_mol13144_rep2_completed.tar.gz \
  md_runs/production/sglt2_mol13144/rep2
sha256sum /workspace/sglt2_mol13144_rep2_completed.tar.gz \
  > /workspace/sglt2_mol13144_rep2_completed.tar.gz.sha256
```

Download and verify that archive before eventually terminating the pod.

## Fetch and extract the Mol_13733 package

```bash
cd /workspace/anpdb-workflow
git pull --ff-only origin agent/runpod-sglt2-rep2

sha256sum -c \
  transfer/runpod/sglt2_mol13733_reps1-3_portable_*.tar.gz.sha256
tar -xzf transfer/runpod/sglt2_mol13733_reps1-3_portable_*.tar.gz
```

The existing `/workspace/miniforge3/envs/anpdb-md` environment is reused.

## Preflight all three replicates

```bash
cd /workspace/anpdb-workflow
PYTHON=/workspace/miniforge3/envs/anpdb-md/bin/python \
  scripts/run_runpod_sglt2_mol13733.sh --check-only
```

The output must report `CUDA context: OK`, `input_ready: true` for every
config, and `RunPod preflight passed for SGLT2 Mol_13733 replicates 1-3`.

## Run sequentially in tmux

```bash
tmux new-session -s sglt2_mol13733
cd /workspace/anpdb-workflow
PYTHON=/workspace/miniforge3/envs/anpdb-md/bin/python \
  scripts/run_runpod_sglt2_mol13733.sh
```

Detach with `Ctrl-b`, then `d`. The launcher automatically starts replicate 2
after replicate 1 completes, and replicate 3 after replicate 2 completes.

Monitor the active replicate with:

```bash
tail -n 3 md_runs/production/sglt2_mol13733/rep*/production.log
nvidia-smi
```

## Recovery after a pod or SSH interruption

If only SSH disconnected, reattach with:

```bash
tmux attach -t sglt2_mol13733
```

If the pod process stopped, start a new tmux session and rerun the same
launcher. It skips completed replicates and uses `--resume` for a partial
replicate when both production output and a RunPod checkpoint exist. Do not
delete or replace a partial replicate directory.

Checkpoints are written every 500 ps. If a pod stops before the first RunPod
checkpoint is written, contact the workflow maintainer before modifying the
partial files.

## Final retrieval

After all three complete, verify that every production log ends at step
`53750000`, then archive the whole compound directory:

```bash
tar -czf /workspace/sglt2_mol13733_reps1-3_completed.tar.gz \
  md_runs/production/sglt2_mol13733
sha256sum /workspace/sglt2_mol13733_reps1-3_completed.tar.gz \
  > /workspace/sglt2_mol13733_reps1-3_completed.tar.gz.sha256
```

Download both files and verify the checksum before terminating the pod.
