#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/run_pending_md_equilibrations.sh
  scripts/run_pending_md_equilibrations.sh --dry-run-only
  scripts/run_pending_md_equilibrations.sh --force
  scripts/run_pending_md_equilibrations.sh --archive

Runs the pending corrected primary MD equilibrations:
  - SGLT2 / Mol_13144, Mol_13733, Mol_15088
  - OPRK1 / Mol_16614
  - 3 replicates each

The script is suitable for RunPod/Linux and WSL. On WSL it uses
scripts/run_openmm_wsl.sh automatically; otherwise it uses ${PYTHON:-python}.

Options:
  --dry-run-only   Check all configs but do not run minimisation/equilibration.
  --force          Run even if equilibration_manifest.json and production.chk exist.
  --no-dry-run     Skip the preflight dry-run before each equilibration.
  --archive        Create a review tarball after finishing.
  --python PATH    Python executable to use on non-WSL Linux.
  -h, --help       Show this message.
EOF
}

dry_run_only=false
force=false
preflight=true
archive=false
python_bin="${PYTHON:-python}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run-only)
      dry_run_only=true
      ;;
    --force)
      force=true
      ;;
    --no-dry-run)
      preflight=false
      ;;
    --archive)
      archive=true
      ;;
    --python)
      if [[ $# -lt 2 ]]; then
        echo "--python requires a path" >&2
        exit 2
      fi
      python_bin="$2"
      shift
      ;;
    -h|--help|help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
  shift
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

configs=(
  configs/md/production/sglt2_mol13144_100ns_rep1.toml
  configs/md/production/sglt2_mol13144_100ns_rep2.toml
  configs/md/production/sglt2_mol13144_100ns_rep3.toml
  configs/md/production/sglt2_mol13733_100ns_rep1.toml
  configs/md/production/sglt2_mol13733_100ns_rep2.toml
  configs/md/production/sglt2_mol13733_100ns_rep3.toml
  configs/md/production/sglt2_mol15088_100ns_rep1.toml
  configs/md/production/sglt2_mol15088_100ns_rep2.toml
  configs/md/production/sglt2_mol15088_100ns_rep3.toml
  configs/md/production/oprk1_mol16614_100ns_rep1.toml
  configs/md/production/oprk1_mol16614_100ns_rep2.toml
  configs/md/production/oprk1_mol16614_100ns_rep3.toml
)

if grep -qi microsoft /proc/version 2>/dev/null; then
  runner=(scripts/run_openmm_wsl.sh)
elif "${python_bin}" -c 'import tomllib' >/dev/null 2>&1; then
  runner=("${python_bin}")
elif command -v python3 >/dev/null 2>&1 && python3 -c 'import tomllib' >/dev/null 2>&1; then
  runner=(python3)
else
  cat >&2 <<'EOF'
No suitable Python was found. Activate the MD environment first, or pass:
  --python /path/to/env/bin/python

The selected Python must be 3.11+ or otherwise provide tomllib.
EOF
  exit 1
fi

log_dir="md_runs/equilibration_logs"
mkdir -p "${log_dir}"
timestamp="$(date -u +%Y%m%d_%H%M%S)"

if [[ "${dry_run_only}" == false ]]; then
  echo "== Checking MD Python environment"
  "${runner[@]}" -c 'import openmm, pdbfixer; import openff.toolkit, openmmforcefields; print("MD imports OK")'
  if [[ -f scripts/check_openmm_cuda.py ]]; then
    echo "== Checking OpenMM CUDA platform"
    "${runner[@]}" scripts/check_openmm_cuda.py --device-index 0 --precision mixed
  fi
fi

config_output_dir() {
  local config="$1"
  sed -nE 's/^[[:space:]]*output_dir[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p' "${config}" | head -n 1
}

run_config() {
  local config="$1"
  local output_dir
  output_dir="$(config_output_dir "${config}")"
  if [[ -z "${output_dir}" ]]; then
    echo "Could not find run.output_dir in ${config}" >&2
    exit 1
  fi

  local run_dir="${repo_root}/${output_dir}"
  local run_name
  run_name="$(basename "$(dirname "${run_dir}")")_$(basename "${run_dir}")"
  local log_path="${log_dir}/${timestamp}_${run_name}_equilibrate.log"

  if [[ "${force}" == false && -f "${run_dir}/equilibration_manifest.json" && -f "${run_dir}/production.chk" ]]; then
    echo "SKIP ${config}: equilibration already present at ${run_dir}"
    return
  fi

  if [[ "${preflight}" == true ]]; then
    echo
    echo "== Dry-run: ${config}"
    "${runner[@]}" scripts/md_production.py --config "${config}" --dry-run
  fi

  if [[ "${dry_run_only}" == true ]]; then
    return
  fi

  echo
  echo "== Equilibrating: ${config}"
  echo "== Log: ${log_path}"
  "${runner[@]}" scripts/md_production.py --config "${config}" --equilibrate-only 2>&1 | tee "${log_path}"
}

for config in "${configs[@]}"; do
  run_config "${config}"
done

if [[ "${archive}" == true && "${dry_run_only}" == false ]]; then
  archive_path="md_equilibration_review_${timestamp}.tar.gz"
  tar -czf "${archive_path}" \
    md_runs/equilibration_logs \
    md_runs/production/sglt2_mol13144 \
    md_runs/production/sglt2_mol13733 \
    md_runs/production/sglt2_mol15088 \
    md_runs/production/oprk1_mol16614
  echo
  echo "Wrote ${archive_path}"
fi

echo
echo "Done."
