#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/ElectricalResistivity"

if [[ ! -x "${BIN}" ]]; then
  echo "Error: binary not found or not executable: ${BIN}" >&2
  exit 1
fi

LIB_PATHS="/home/bcellier/Documents/Stage/new_install/pagure/metis/gcc133/5.1.0/lib:/home/bcellier/Documents/Stage/new_install/pagure/suitesparse/gcc133/5.10.1/lib:/home/bcellier/Documents/Stage/new_install/pagure/cgal/gcc133/5.5.2/lib:/home/bcellier/Documents/Stage/new_install/pagure/boost/py37/gcc133/1.84.0/lib"

# Keep existing LD_LIBRARY_PATH if already defined.
if [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
  export LD_LIBRARY_PATH="${LIB_PATHS}:${LD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="${LIB_PATHS}"
fi

exec "${BIN}" "$@"
