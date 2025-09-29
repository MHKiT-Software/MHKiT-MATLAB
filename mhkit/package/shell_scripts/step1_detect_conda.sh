#!/bin/bash

# MHKiT-MATLAB Step 1: Conda Detection
# Generated from template - DO NOT EDIT DIRECTLY
# Report issues: https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/new

set -e  # Exit on any error

# Script configuration
CONDA_ENV_NAME="mhkit-matlab-env"
PYTHON_VERSION="3.11"
MHKIT_VERSION="0.9"
GITHUB_ISSUES="https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/new"

# Logging functions
log_info() {
    echo "[INFO] $1"
}

log_success() {
    echo "[SUCCESS] $1"
}

log_warning() {
    echo "[WARNING] $1"
}

log_error() {
    echo "[ERROR] $1"
    echo "[ERROR] Please report this issue: $GITHUB_ISSUES"
}

# Function to expand path with environment variables and ~
expand_path() {
    local path="$1"
    # Expand ~ to $HOME
    path="${path/#\~/$HOME}"
    # Expand environment variables
    path=$(eval echo "$path")
    echo "$path"
}

# Function to test if a conda installation works
test_conda_installation() {
    local conda_path="$1"
    log_info "Testing conda at: $conda_path"

    # Test if the conda executable exists and is executable
    if [[ -x "$conda_path" ]]; then
        # Test if conda info works (this is the most reliable test)
        if "$conda_path" info >/dev/null 2>&1; then
            local conda_version=$("$conda_path" --version 2>/dev/null || echo "unknown")
            local conda_location=$("$conda_path" info --base 2>/dev/null || echo "unknown")
            log_success "Working conda found!"
            log_info "  Version: $conda_version"
            log_info "  Location: $conda_location"
            return 0
        else
            log_warning "Conda executable found but 'conda info' failed"
            return 1
        fi
    else
        return 1
    fi
}

# OS Detection
OS_TYPE=$(uname)
ARCH_TYPE=$(uname -m)
# Function to check if conda is available
detect_conda() {
    log_info "Starting comprehensive conda detection..."

    # Step 1: Check if conda is already in PATH and working
    if command -v conda >/dev/null 2>&1; then
        log_info "Found conda in PATH"
        if test_conda_installation "$(command -v conda)"; then
            echo "CONDA_DETECTED=true"
            echo "CONDA_PATH=$(command -v conda)"
            return 0
        else
            log_warning "Conda in PATH is not working properly"
        fi
    fi

    # Step 2: Check environment variables
    if [[ -n "$${CONDA_EXE:-}" ]]; then
        log_info "Checking environment variable CONDA_EXE: $${CONDA_EXE}"
        local conda_from_env=""
        if [[ "CONDA_EXE" == "CONDA_EXE" ]]; then
            conda_from_env="$${CONDA_EXE}"
        elif [[ "CONDA_EXE" == "CONDA_PREFIX" ]]; then
            conda_from_env="$${CONDA_EXE}/bin/conda"
        fi

        if [[ -n "$conda_from_env" ]]; then
            conda_from_env=$(expand_path "$conda_from_env")
            if test_conda_installation "$conda_from_env"; then
                echo "CONDA_DETECTED=true"
                echo "CONDA_PATH=$conda_from_env"
                return 0
            fi
        fi
    fi
    if [[ -n "$${CONDA_PREFIX:-}" ]]; then
        log_info "Checking environment variable CONDA_PREFIX: $${CONDA_PREFIX}"
        local conda_from_env=""
        if [[ "CONDA_PREFIX" == "CONDA_EXE" ]]; then
            conda_from_env="$${CONDA_PREFIX}"
        elif [[ "CONDA_PREFIX" == "CONDA_PREFIX" ]]; then
            conda_from_env="$${CONDA_PREFIX}/bin/conda"
        fi

        if [[ -n "$conda_from_env" ]]; then
            conda_from_env=$(expand_path "$conda_from_env")
            if test_conda_installation "$conda_from_env"; then
                echo "CONDA_DETECTED=true"
                echo "CONDA_PATH=$conda_from_env"
                return 0
            fi
        fi
    fi

    # Step 3: Check all common installation paths systematically
    log_info "Checking common conda installation paths..."
    local paths_to_check=(
        "~/miniconda3/bin/conda"
        "~/miniconda3/condabin/conda"
        "~/anaconda3/bin/conda"
        "~/anaconda3/condabin/conda"
        "/opt/miniconda3/bin/conda"
        "/opt/miniconda3/condabin/conda"
        "/opt/anaconda3/bin/conda"
        "/opt/anaconda3/condabin/conda"
        "/usr/local/miniconda3/bin/conda"
        "/usr/local/miniconda3/condabin/conda"
        "/usr/local/anaconda3/bin/conda"
        "/usr/local/anaconda3/condabin/conda"
        "/usr/bin/conda"
        "/usr/local/bin/conda"
        "~/.local/bin/conda"
        "/opt/homebrew/bin/conda"
        "/opt/homebrew/condabin/conda"
        "/usr/local/Caskroom/miniconda/base/bin/conda"
        "/usr/local/Caskroom/miniconda/base/condabin/conda"
    )

    for path in "${paths_to_check[@]}"; do
        local expanded_path=$(expand_path "$path")
        if test_conda_installation "$expanded_path"; then
            echo "CONDA_DETECTED=true"
            echo "CONDA_PATH=$expanded_path"
            return 0
        fi
    done

    log_info "No working conda installation found"
    echo "CONDA_DETECTED=false"
    return 1
}

# Main execution
main() {
    log_info "Step 1: Detecting conda installation..."
    detect_conda
}

# Run main function
main "$@"