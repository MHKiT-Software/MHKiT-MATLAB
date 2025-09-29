# MHKiT-MATLAB Step 5: Post-Install Configuration
# Generated from template - DO NOT EDIT DIRECTLY
# Report issues: https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/new

$ErrorActionPreference = "Stop"  # Exit on any error

# Script configuration
$CONDA_ENV_NAME = "mhkit-matlab-env"
$PYTHON_VERSION = "3.11"
$MHKIT_VERSION = "0.9"
$GITHUB_ISSUES = "https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/new"

# Logging functions
function Write-Info {
    param([string]$Message)
    Write-Host "[INFO] $Message" -ForegroundColor Green
}

function Write-Success {
    param([string]$Message)
    Write-Host "[SUCCESS] $Message" -ForegroundColor Cyan
}

function Write-Warning {
    param([string]$Message)
    Write-Host "[WARNING] $Message" -ForegroundColor Yellow
}

function Write-Error {
    param([string]$Message)
    Write-Host "[ERROR] $Message" -ForegroundColor Red
    Write-Host "[ERROR] Please report this issue: $GITHUB_ISSUES" -ForegroundColor Red
}

# Function to test if a conda installation works
function Test-Conda {
    param([string]$CondaPath)

    Write-Info "Testing conda at: $CondaPath"

    # Test if the conda executable exists
    if (Test-Path $CondaPath) {
        try {
            # Test if conda info works (this is the most reliable test)
            $null = & "$CondaPath" info 2>$null
            $condaVersion = & "$CondaPath" --version 2>$null
            $condaLocation = & "$CondaPath" info --base 2>$null
            Write-Success "Working conda found!"
            Write-Info "  Version: $condaVersion"
            Write-Info "  Location: $condaLocation"
            return $true
        } catch {
            Write-Warning "Conda executable found but 'conda info' failed"
            return $false
        }
    } else {
        return $false
    }
}

# OS and Architecture Detection
$OSType = "Windows"
$ArchType = $env:PROCESSOR_ARCHITECTURE
# Function to install post-install dependencies
function Install-PostDependencies {
    param([string]$CondaPath)

    Write-Info "Installing post-install dependencies..."

    Write-Info "Executing: conda run -n mhkit-matlab-env pip install --upgrade netcdf4"
    try {
        conda run -n mhkit-matlab-env pip install --upgrade netcdf4
    } catch {
        Write-Error "Failed to execute: conda run -n mhkit-matlab-env pip install --upgrade netcdf4"
        Write-Output "POST_DEPS_INSTALLED=false"
        exit 1
    }
    Write-Info "Executing: conda run -n mhkit-matlab-env conda install scipy -y"
    try {
        conda run -n mhkit-matlab-env conda install scipy -y
    } catch {
        Write-Error "Failed to execute: conda run -n mhkit-matlab-env conda install scipy -y"
        Write-Output "POST_DEPS_INSTALLED=false"
        exit 1
    }

    Write-Success "Post-install dependencies installed"
    Write-Output "POST_DEPS_INSTALLED=true"
}

# Function to test functionality
function Test-Functionality {
    param([string]$CondaPath)

    Write-Info "Testing mhkit functionality..."

    try {
        $testResult = & "$CondaPath" run -n $CONDA_ENV_NAME python -c "import mhkit; print(mhkit.river.performance.circular(30))" 2>$null

        if ($testResult -like "*30*" -and $testResult -like "*706*") {
            Write-Success "Functionality test passed: $testResult"
            Write-Output "FUNCTIONALITY_TEST=passed"
        } else {
            Write-Error "Functionality test failed. Output: $testResult"
            Write-Output "FUNCTIONALITY_TEST=failed"
            exit 1
        }
    } catch {
        Write-Error "Functionality test failed: $_"
        Write-Output "FUNCTIONALITY_TEST=failed"
        exit 1
    }
}

# Function to get Python executable path
function Get-PythonPath {
    param([string]$CondaPath)

    try {
        $pythonPath = & "$CondaPath" run -n $CONDA_ENV_NAME python -c "import sys; print(sys.executable)"
        Write-Output "PYTHON_PATH=$pythonPath"
    } catch {
        Write-Error "Failed to get Python path: $_"
        exit 1
    }
}

# Main execution
function Main {
    param([string]$CondaPath)

    if (-not $CondaPath) {
        Write-Error "No conda path provided. Usage: script.ps1 <conda_path>"
        exit 1
    }

    Write-Info "Step 5: Post-install configuration..."
    Write-Info "Using conda at: $CondaPath"

    Install-PostDependencies $CondaPath
    Test-Functionality $CondaPath
    Get-PythonPath $CondaPath

    Write-Success "Post-install configuration completed!"
}

# Run main function with arguments
Main $args[0]