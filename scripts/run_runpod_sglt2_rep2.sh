#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

python_bin="${PYTHON:-python}"
config="configs/md/production/sglt2_mol13144_100ns_rep2.toml"
run_dir="md_runs/production/sglt2_mol13144/rep2"
check_only=false
[[ "${1:-}" == "--check-only" ]] && check_only=true

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

if [[ -s "${run_dir}/production.dcd" || -s "${run_dir}/production.log" ]]; then
  echo "Refusing portable restart because production output already exists in ${run_dir}" >&2
  exit 1
fi

"${python_bin}" -c 'import openmm, pdbfixer; import openff.toolkit, openmmforcefields; print("MD imports OK; OpenMM", openmm.version.version)'
"${python_bin}" scripts/check_openmm_cuda.py --device-index 0 --precision mixed
"${python_bin}" scripts/md_production.py --config "${config}" --dry-run

if [[ "${check_only}" == true ]]; then
  echo "RunPod preflight passed."
  exit 0
fi

mkdir -p "${run_dir}"
echo "Starting SGLT2 Mol_13144 replicate 2 from portable equilibrated_state.xml"
"${python_bin}" scripts/md_production.py --config "${config}" --resume-state \
  2>&1 | tee "${run_dir}/runpod_production_console.log"

