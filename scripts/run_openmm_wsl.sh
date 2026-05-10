#!/usr/bin/env bash
set -euo pipefail

REAL_HOME="$(getent passwd "${USER}" | cut -d: -f6)"
if [[ -z "${REAL_HOME}" ]]; then
  REAL_HOME="${HOME}"
fi

CONDA_ROOT="${CONDA_ROOT:-$REAL_HOME/miniforge3}"
ENV_NAME="${OPENMM_ENV_NAME:-anpdb-md}"
ENV_PREFIX="${CONDA_ROOT}/envs/${ENV_NAME}"

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  echo "OpenMM environment not found at ${ENV_PREFIX}" >&2
  exit 1
fi

if [[ -n "${SNAP_NAME:-}" || "${PATH:-}" == *"/snap/"* ]]; then
  cat >&2 <<'EOF'
This launcher must be run from a normal WSL shell, not from inside the codex snap.
Open a regular Ubuntu terminal and rerun the same command there.
EOF
  exit 2
fi

if [[ -d /usr/lib/wsl/lib ]]; then
  WSL_GPU_LIB="/usr/lib/wsl/lib"
elif [[ -d /var/lib/snapd/hostfs/usr/lib/wsl/lib ]]; then
  WSL_GPU_LIB="/var/lib/snapd/hostfs/usr/lib/wsl/lib"
else
  echo "WSL GPU driver library directory not found" >&2
  exit 1
fi

export OPENMM_PLUGIN_DIR="${ENV_PREFIX}/lib/plugins"
export LD_LIBRARY_PATH="${WSL_GPU_LIB}:${ENV_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

exec "${ENV_PREFIX}/bin/python" "$@"
