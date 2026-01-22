#!/bin/bash

# MHKiT-MATLAB Step 3: Create Conda Environment
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
# Function to check if conda environment exists
check_conda_env() {
    local conda_path="$1"
    if "$conda_path" env list | grep -q "^$CONDA_ENV_NAME "; then
        log_success "Found existing conda environment: $CONDA_ENV_NAME"
        echo "ENV_EXISTS=true"
        return 0
    else
        log_info "Conda environment '$CONDA_ENV_NAME' not found, will create it"
        echo "ENV_EXISTS=false"
        return 1
    fi
}

# Function to create conda environment
create_conda_env() {
    local conda_path="$1"
    log_info "Creating conda environment '$CONDA_ENV_NAME' with Python $PYTHON_VERSION..."

    # Accept conda Terms of Service to avoid interactive prompts
    log_info "Accepting conda Terms of Service..."
    "$conda_path" config --set solver libmamba 2>/dev/null || true

    # Try to accept TOS (may fail on older conda versions, that's OK)
    "$conda_path" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    "$conda_path" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

    # Create environment using conda-forge to minimize TOS issues
    "$conda_path" create -n "$CONDA_ENV_NAME" python="$PYTHON_VERSION" -c conda-forge -y

    if [ $? -eq 0 ]; then
        log_success "Conda environment created successfully"
        echo "ENV_CREATED=true"
    else
        log_error "Failed to create conda environment"
        echo "ENV_CREATED=false"
        exit 1
    fi
}

# Function to verify Python version in environment
verify_python_version() {
    local conda_path="$1"
    log_info "Verifying Python version in environment..."

    local py_version=$("$conda_path" run -n "$CONDA_ENV_NAME" python --version 2>&1 | cut -d' ' -f2)
    local major_minor=$(echo "$py_version" | cut -d'.' -f1,2)
    local expected_major_minor=$(echo "$PYTHON_VERSION" | cut -d'.' -f1,2)

    if [[ "$major_minor" == "$expected_major_minor" ]]; then
        log_success "Python version verified: $py_version"
        echo "PYTHON_VERSION_OK=true"
        return 0
    else
        log_warning "Python version mismatch. Expected: $PYTHON_VERSION, Got: $py_version"
        log_info "Environment will need to be recreated..."
        echo "PYTHON_VERSION_OK=false"
        return 1
    fi
}

main() {
    local conda_path="$1"
    if [[ -z "$conda_path" ]]; then
        log_error "No conda path provided. Usage: $0 <conda_path>"
        exit 1
    fi

    log_info "Step 3: Managing conda environment..."
    log_info "Using conda at: $conda_path"

    # Add conda to PATH for this session
    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    if ! check_conda_env "$conda_path"; then
        create_conda_env "$conda_path"
    fi

    if ! verify_python_version "$conda_path"; then
        log_info "Recreating environment with correct Python version..."
        "$conda_path" remove -n "$CONDA_ENV_NAME" --all -y >/dev/null 2>&1 || true
        create_conda_env "$conda_path"
        verify_python_version "$conda_path"
    fi
}

main "$@"
