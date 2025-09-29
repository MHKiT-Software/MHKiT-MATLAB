#!/bin/bash

# MHKiT-MATLAB Step 4: Install Dependencies
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
# Function to install pre-install dependencies
install_pre_dependencies() {
    local conda_path="$1"
    log_info "Installing pre-install dependencies..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    if [[ "$OS_TYPE" == "Darwin" ]]; then
        log_info "Executing: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y"
        conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y"
            echo "PRE_DEPS_INSTALLED=false"
            exit 1
        fi
    else
        log_info "Executing: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y"
        conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf cftime netcdf4 pandas numpy -c conda-forge -y"
            echo "PRE_DEPS_INSTALLED=false"
            exit 1
        fi
    fi

    log_success "Pre-install dependencies installed"
    echo "PRE_DEPS_INSTALLED=true"
}

# Function to check if mhkit is already installed with correct version
check_mhkit_installed() {
    local conda_path="$1"
    log_info "Checking if mhkit-python is already installed..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    local installed_version=$("$conda_path" run -n "$CONDA_ENV_NAME" python -c "import mhkit; print(mhkit.__version__)" 2>/dev/null || echo "not_installed")

    if [[ "$installed_version" == *"$MHKIT_VERSION"* ]]; then
        log_success "mhkit-python v$MHKIT_VERSION already installed"
        echo "MHKIT_INSTALLED=true"
        return 0
    else
        log_info "mhkit-python not installed or wrong version (found: $installed_version)"
        echo "MHKIT_INSTALLED=false"
        return 1
    fi
}

# Function to install mhkit-python
install_mhkit() {
    local conda_path="$1"
    log_info "Installing mhkit-python v$MHKIT_VERSION..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    "$conda_path" run -n "$CONDA_ENV_NAME" conda install -c conda-forge mhkit=="$MHKIT_VERSION" -y

    if [ $? -eq 0 ]; then
        log_success "mhkit-python installed successfully"
        echo "MHKIT_INSTALLED=true"
    else
        log_error "Failed to install mhkit-python"
        echo "MHKIT_INSTALLED=false"
        exit 1
    fi

    # Verify installation
    log_info "Verifying mhkit-python installation..."
    local installed_version=$("$conda_path" run -n "$CONDA_ENV_NAME" python -c "import mhkit; print(mhkit.__version__)")

    if [[ "$installed_version" == *"$MHKIT_VERSION"* ]]; then
        log_success "Installation verified: $installed_version"
        echo "MHKIT_VERIFIED=true"
    else
        log_error "Version verification failed. Expected: v$MHKIT_VERSION, Got: $installed_version"
        echo "MHKIT_VERIFIED=false"
        exit 1
    fi
}

# Function to install mhkit_python_utils
install_mhkit_python_utils() {
    local conda_path="$1"
    log_info "Installing mhkit_python_utils..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    # Create temporary directory for download
    local temp_dir=$(mktemp -d)
    local utils_url="https://github.com/MHKiT-Software/MHKiT-MATLAB/releases/download/v0.6.0/mhkit_python_utils_v0.6.0.zip"

    log_info "Downloading mhkit_python_utils from: $utils_url"

    # Download the zip file
    if command -v curl >/dev/null 2>&1; then
        curl -L "$utils_url" -o "$temp_dir/mhkit_python_utils.zip"
    elif command -v wget >/dev/null 2>&1; then
        wget "$utils_url" -O "$temp_dir/mhkit_python_utils.zip"
    else
        log_error "Neither curl nor wget found. Cannot download mhkit_python_utils."
        echo "UTILS_INSTALLED=false"
        exit 1
    fi

    # Unzip the package
    log_info "Extracting mhkit_python_utils..."
    cd "$temp_dir"
    unzip -q mhkit_python_utils.zip

    # Find the extracted directory (should be mhkit_python_utils or similar)
    local utils_dir=$(find . -name "*mhkit_python_utils*" -type d | head -1)
    if [[ -z "$utils_dir" ]]; then
        log_error "Could not find extracted mhkit_python_utils directory"
        echo "UTILS_INSTALLED=false"
        exit 1
    fi

    # Install in development mode
    log_info "Installing mhkit_python_utils with pip install -e ..."
    cd "$utils_dir"
    "$conda_path" run -n "$CONDA_ENV_NAME" pip install -e .

    if [ $? -eq 0 ]; then
        log_success "mhkit_python_utils installed successfully"
        echo "UTILS_INSTALLED=true"

        # Verify installation
        log_info "Verifying mhkit_python_utils installation..."
        local utils_version=$("$conda_path" run -n "$CONDA_ENV_NAME" python -c "import mhkit_python_utils; print(mhkit_python_utils.__version__)" 2>/dev/null || echo "unknown")
        log_info "mhkit_python_utils version: $utils_version"

        # Clean up temporary directory
        rm -rf "$temp_dir"
    else
        log_error "Failed to install mhkit_python_utils"
        echo "UTILS_INSTALLED=false"
        rm -rf "$temp_dir"
        exit 1
    fi
}

main() {
    local conda_path="$1"
    if [[ -z "$conda_path" ]]; then
        log_error "No conda path provided. Usage: $0 <conda_path>"
        exit 1
    fi

    log_info "Step 4: Installing dependencies..."
    log_info "Using conda at: $conda_path"

    install_pre_dependencies "$conda_path"

    if ! check_mhkit_installed "$conda_path"; then
        install_mhkit "$conda_path"
    fi

    install_mhkit_python_utils "$conda_path"
}

main "$@"
