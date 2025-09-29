#!/bin/bash

# MHKiT-MATLAB Step 5: Post-Install Configuration
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
# Function to install post-install dependencies
install_post_dependencies() {
    local conda_path="$1"
    log_info "Installing post-install dependencies..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    if [[ "$OS_TYPE" == "Darwin" ]]; then
        log_info "Executing: conda run -n mhkit-matlab-env conda install openssl=3.0.* ca-certificates -c conda-forge -y"
        conda run -n mhkit-matlab-env conda install openssl=3.0.* ca-certificates -c conda-forge -y
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env conda install openssl=3.0.* ca-certificates -c conda-forge -y"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env pip uninstall -y netcdf4"
        conda run -n mhkit-matlab-env pip uninstall -y netcdf4
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip uninstall -y netcdf4"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env pip install netcdf4 --no-deps --no-cache-dir"
        conda run -n mhkit-matlab-env pip install netcdf4 --no-deps --no-cache-dir
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip install netcdf4 --no-deps --no-cache-dir"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env pip install cftime --no-cache-dir"
        conda run -n mhkit-matlab-env pip install cftime --no-cache-dir
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip install cftime --no-cache-dir"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env pip uninstall -y scipy"
        conda run -n mhkit-matlab-env pip uninstall -y scipy
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip uninstall -y scipy"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env pip install scipy==1.12.0 --no-cache-dir"
        conda run -n mhkit-matlab-env pip install scipy==1.12.0 --no-cache-dir
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip install scipy==1.12.0 --no-cache-dir"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env conda clean --all -y"
        conda run -n mhkit-matlab-env conda clean --all -y
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env conda clean --all -y"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
    else
        log_info "Executing: conda run -n mhkit-matlab-env pip install --upgrade netcdf4"
        conda run -n mhkit-matlab-env pip install --upgrade netcdf4
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env pip install --upgrade netcdf4"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
        log_info "Executing: conda run -n mhkit-matlab-env conda install scipy -y"
        conda run -n mhkit-matlab-env conda install scipy -y
        if [ $? -ne 0 ]; then
            log_error "Failed to execute: conda run -n mhkit-matlab-env conda install scipy -y"
            echo "POST_DEPS_INSTALLED=false"
            exit 1
        fi
    fi

    log_success "Post-install dependencies installed"
    echo "POST_DEPS_INSTALLED=true"
}

# Function to test functionality
test_functionality() {
    local conda_path="$1"
    log_info "Testing mhkit functionality..."

    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    local test_result=$("$conda_path" run -n "$CONDA_ENV_NAME" python -c "import mhkit; print(mhkit.river.performance.circular(30))" 2>/dev/null)

    if [[ "$test_result" == *"30"* ]] && [[ "$test_result" == *"706"* ]]; then
        log_success "Functionality test passed: $test_result"
        echo "FUNCTIONALITY_TEST=passed"
    else
        log_error "Functionality test failed. Output: $test_result"
        echo "FUNCTIONALITY_TEST=failed"
        exit 1
    fi
}

# Function to get Python executable path
get_python_path() {
    local conda_path="$1"
    local conda_dir=$(dirname "$conda_path")
    export PATH="$conda_dir:$PATH"

    local python_path=$("$conda_path" run -n "$CONDA_ENV_NAME" which python)
    echo "PYTHON_PATH=$python_path"
}

# Main execution
main() {
    local conda_path="$1"
    if [[ -z "$conda_path" ]]; then
        log_error "No conda path provided. Usage: $0 <conda_path>"
        exit 1
    fi

    log_info "Step 5: Post-install configuration..."
    log_info "Using conda at: $conda_path"

    install_post_dependencies "$conda_path"
    test_functionality "$conda_path"
    get_python_path "$conda_path"

    log_success "Post-install configuration completed!"
}

# Run main function
main "$@"