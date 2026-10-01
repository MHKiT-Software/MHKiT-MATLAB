#!/usr/bin/env bash
#
# Run the MATLAB-only test suites locally on macOS
#
#   scripts/run_matlab_only_tests_macos.sh            # unit tests (default)
#   scripts/run_matlab_only_tests_macos.sh tests      # unit tests
#   scripts/run_matlab_only_tests_macos.sh examples   # example live scripts
#   scripts/run_matlab_only_tests_macos.sh all        # both
#
# MATLAB is located from (in order): $MATLAB_EXE, `matlab` on PATH, the newest
# /Applications/MATLAB_R*.app on macOS.
#
# To isolate a users python config MATLAB is started with a minimal PATH
# and a throw-away preferences directory. This means any user `pyenv` MATLAB 
# configuration is ignored

set -euo pipefail

what="${1:-tests}"
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

case "$what" in
    tests)    cmd="results = runMatlabOnlyTests();" ;;
    examples) cmd="results = runMatlabOnlyExamples();" ;;
    all)      cmd="results = [runMatlabOnlyTests(), runMatlabOnlyExamples()];" ;;
    *) echo "Usage: $0 [tests|examples|all]" >&2; exit 2 ;;
esac

matlab_exe="$(find_matlab)"
prefdir="$(mktemp -d "${TMPDIR:-/tmp}/mhkit_matlab_prefs.XXXXXX")"
trap 'rm -rf "$prefdir"' EXIT

echo "MATLAB:   $matlab_exe"
echo "Running:  $what"
echo "Prefdir:  $prefdir (temporary)"

cd "$repo_root"
env -i \
    HOME="$HOME" \
    PATH="/usr/bin:/bin:/usr/sbin:/sbin" \
    MATLAB_PREFDIR="$prefdir" \
    "$matlab_exe" -batch \
    "addpath(fullfile('mhkit','tests')); $cmd disp(table(results)); assertSuccess(results);"
