#!/bin/bash

# MHKiT-MATLAB Step 2: Install Miniconda
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
# Function to install miniconda
install_conda() {
    log_info "Installing miniconda..."

    if [[ "$OS_TYPE" == "Darwin" ]]; then
        if [[ "$ARCH_TYPE" == "arm64" ]]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
        else
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
        fi
    elif [[ "$OS_TYPE" == "Linux" ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    else
        log_error "Unsupported operating system: $OS_TYPE"
        exit 1
    fi

    # Create miniconda directory
    mkdir -p ~/miniconda3

    # Download miniconda installer
    log_info "Downloading miniconda installer..."
    if command -v curl >/dev/null 2>&1; then
        curl -L "$MINICONDA_URL" -o ~/miniconda3/miniconda.sh
    elif command -v wget >/dev/null 2>&1; then
        wget "$MINICONDA_URL" -O ~/miniconda3/miniconda.sh
    else
        log_error "Neither curl nor wget found. Please install one of them first."
        exit 1
    fi

    # Install miniconda
    log_info "Installing miniconda to ~/miniconda3..."
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

    # Remove installer
    rm ~/miniconda3/miniconda.sh

    # Initialize conda for future shells
    ~/miniconda3/bin/conda init bash >/dev/null 2>&1 || true

    log_success "Miniconda installed successfully"

    # Output conda path for next steps
    echo "CONDA_PATH=$HOME/miniconda3/bin/conda"
}

# Main execution
main() {
    log_info "Step 2: Installing miniconda..."
    install_conda
}

# Run main function
main "$@"