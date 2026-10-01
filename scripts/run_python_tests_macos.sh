#!/usr/bin/env bash
#
# Run the full MHKiT-MATLAB unit test suite, including tests that call
# MHKiT-Python, locally on macOS
#
#   scripts/run_python_tests_macos.sh /path/to/python
#
# The Python interpreter must have mhkit and mhkit_python_utils installed:
#
#   conda create -n mhkit_matlab -c conda-forge python=3.12 "numpy>=2" pip netcdf4 hdf5
#   conda activate mhkit_matlab
#   pip install "mhkit[all]==1.1.2" "pandas<3"  # pecos 1.0.0 check_delta fails with pandas 3
#   pip install -e .
#   scripts/run_python_tests_macos.sh "$(which python)"
#
# MATLAB is located from (in order): $MATLAB_EXE, `matlab` on PATH, the newest
# /Applications/MATLAB_R*.app on macOS.
#
# To isolate a users python config MATLAB is started with a minimal PATH
# and a throw-away preferences directory. This means any user `pyenv` MATLAB
# configuration is ignored

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/python" >&2
    exit 2
fi
python_exe="$1"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

find_matlab() {
    if [[ -n "${MATLAB_EXE:-}" ]]; then
        echo "$MATLAB_EXE"; return
    fi
    if command -v matlab >/dev/null 2>&1; then
        command -v matlab; return
    fi
    local newest
    newest="$(ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n 1 || true)"
    if [[ -n "$newest" ]]; then
        echo "$newest/bin/matlab"; return
    fi
    echo "Could not find MATLAB. Set MATLAB_EXE=/path/to/matlab" >&2
    exit 1
}

"$python_exe" -c "import mhkit, mhkit_python_utils; print('MHKiT-Python', mhkit.__version__)"

matlab_exe="$(find_matlab)"
prefdir="$(mktemp -d "${TMPDIR:-/tmp}/mhkit_matlab_prefs.XXXXXX")"
trap 'rm -rf "$prefdir"' EXIT

echo "MATLAB:   $matlab_exe"
echo "Python:   $python_exe"
echo "Prefdir:  $prefdir (temporary)"

cd "$repo_root"
env -i \
    HOME="$HOME" \
    PATH="$(dirname "$python_exe"):/usr/bin:/bin:/usr/sbin:/sbin" \
    MATLAB_PREFDIR="$prefdir" \
    "$matlab_exe" -batch \
    "addpath(fullfile('mhkit','tests')); results = runPythonTests('$python_exe'); disp(table(results)); assertSuccess(results);"
