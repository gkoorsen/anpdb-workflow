#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

python_bin="${PYTHON:-python}"
check_only=false
[[ "${1:-}" == "--check-only" ]] && check_only=true

"${python_bin}" -c 'import openmm, pdbfixer; import openff.toolkit, openmmforcefields; print("MD imports OK; OpenMM", openmm.version.version)'
"${python_bin}" scripts/check_openmm_cuda.py --device-index 0 --precision mixed

for rep in 1 2 3; do
  config="configs/md/production/sglt2_mol13733_100ns_rep${rep}.toml"
  run_dir="md_runs/production/sglt2_mol13733/rep${rep}"
  required=(
    "${run_dir}/system_solvated.pdb"
    "${run_dir}/system.xml"
    "${run_dir}/integrator.xml"
    "${run_dir}/equilibrated_state.xml"
    "${run_dir}/production_plan.json"
    "${run_dir}/endpoint_analysis_manifest.json"
  )
  for file in "${required[@]}"; do
    [[ -s "${file}" ]] || { echo "Missing or empty: ${file}" >&2; exit 1; }
  done
  "${python_bin}" scripts/md_production.py --config "${config}" --dry-run
done

if [[ "${check_only}" == true ]]; then
  echo "RunPod preflight passed for SGLT2 Mol_13733 replicates 1-3."
  exit 0
fi

for rep in 1 2 3; do
  config="configs/md/production/sglt2_mol13733_100ns_rep${rep}.toml"
  run_dir="md_runs/production/sglt2_mol13733/rep${rep}"
  if [[ -s "${run_dir}/run_manifest.json" && -s "${run_dir}/final_state.xml" ]]; then
    echo "SKIP completed SGLT2 Mol_13733 replicate ${rep}"
    continue
  fi

  if [[ -s "${run_dir}/production.dcd" || -s "${run_dir}/production.log" ]]; then
    [[ -s "${run_dir}/production.chk" ]] || {
      echo "Replicate ${rep} has production output but no checkpoint; refusing unsafe restart" >&2
      exit 1
    }
    mode="--resume"
    echo "Resuming SGLT2 Mol_13733 replicate ${rep} from its RunPod checkpoint"
  else
    mode="--resume-state"
    echo "Starting SGLT2 Mol_13733 replicate ${rep} from portable equilibrated_state.xml"
  fi

  "${python_bin}" scripts/md_production.py --config "${config}" "${mode}" \
    2>&1 | tee "${run_dir}/runpod_production_console.log"
done

echo "Completed all SGLT2 Mol_13733 production replicates"

