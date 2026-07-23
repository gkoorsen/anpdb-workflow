#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

run_dir="md_runs/production/sglt2_mol13144/rep2"
stamp="$(date -u +%Y%m%d_%H%M%S)"
output_dir="transfer/runpod"
archive="${output_dir}/sglt2_mol13144_rep2_portable_${stamp}.tar.gz"

files=(
  environment-md.yml
  configs/md/production/sglt2_mol13144_100ns_rep2.toml
  scripts/check_openmm_cuda.py
  scripts/md_production.py
  scripts/run_runpod_sglt2_rep2.sh
  docs/runpod_sglt2_rep2.md
  data/md_inputs/anpdb_truly_novel_std.csv
  data/md_inputs/poses/Mol_13144_7VSI_opm_oriented_out.pdbqt
  data/md_inputs/receptors/7VSI_opm_oriented_clean.pdb
  "${run_dir}/system_solvated.pdb"
  "${run_dir}/system.xml"
  "${run_dir}/integrator.xml"
  "${run_dir}/equilibrated.pdb"
  "${run_dir}/equilibrated_state.xml"
  "${run_dir}/equilibration_manifest.json"
  "${run_dir}/production_plan.json"
  "${run_dir}/endpoint_analysis_manifest.json"
)

for file in "${files[@]}"; do
  [[ -s "${file}" ]] || { echo "Missing or empty required file: ${file}" >&2; exit 1; }
done

mkdir -p "${output_dir}"
tar -czf "${archive}" "${files[@]}"
sha256sum "${archive}" > "${archive}.sha256"

echo "Created ${archive}"
du -h "${archive}"
echo "Checksum: ${archive}.sha256"
echo "Upload both files to the RunPod persistent volume."

