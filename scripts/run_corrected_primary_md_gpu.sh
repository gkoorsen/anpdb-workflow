#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/run_corrected_primary_md_gpu.sh dry-run
  scripts/run_corrected_primary_md_gpu.sh equilibrate
  scripts/run_corrected_primary_md_gpu.sh resume

Runs the corrected primary OpenFF/OpenMM membrane MD set:
  - SGLT2 / Mol_13144, Mol_13733, Mol_15088
  - OPRK1 / Mol_16614

Run this from a normal WSL shell on the GPU machine. The helper delegates to
scripts/run_openmm_wsl.sh, which uses OPENMM_ENV_NAME=anpdb-md by default.
EOF
}

phase="${1:-dry-run}"
case "${phase}" in
  dry-run)
    md_args=(--dry-run)
    ;;
  equilibrate)
    md_args=(--equilibrate-only)
    ;;
  resume)
    md_args=(--resume)
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

exec scripts/run_openmm_wsl.sh scripts/md_batch.py "${configs[@]}" "${md_args[@]}"
