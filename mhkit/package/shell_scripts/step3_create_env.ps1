# MHKiT-MATLAB Step 3: Create Conda Environment
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
# Function to check if conda environment exists
function Test-CondaEnv {
    param([string]$CondaPath)

    $envList = & "$CondaPath" env list 2>$null
    if ($envList -match "^$CONDA_ENV_NAME ") {
        Write-Success "Found existing conda environment: $CONDA_ENV_NAME"
        Write-Output "ENV_EXISTS=true"
        return $true
    } else {
        Write-Info "Conda environment '$CONDA_ENV_NAME' not found, will create it"
        Write-Output "ENV_EXISTS=false"
        return $false
    }
}

# Function to create conda environment
function New-CondaEnv {
    param([string]$CondaPath)

    Write-Info "Creating conda environment '$CONDA_ENV_NAME' with Python $PYTHON_VERSION..."

    try {
        # Accept conda Terms of Service to avoid interactive prompts
        Write-Info "Accepting conda Terms of Service..."
        & "$CondaPath" config --set solver libmamba 2>$null

        # Create environment using conda-forge to minimize TOS issues
        & "$CondaPath" create -n $CONDA_ENV_NAME python=$PYTHON_VERSION -c conda-forge -y

        Write-Success "Conda environment created successfully"
        Write-Output "ENV_CREATED=true"
    } catch {
        Write-Error "Failed to create conda environment: $_"
        Write-Output "ENV_CREATED=false"
        exit 1
    }
}

# Function to verify Python version in environment
function Test-PythonVersion {
    param([string]$CondaPath)

    Write-Info "Verifying Python version in environment..."

    try {
        $pyVersionOutput = & "$CondaPath" run -n $CONDA_ENV_NAME python --version 2>&1
        $pyVersion = ($pyVersionOutput -split ' ')[1]
        $majorMinor = ($pyVersion -split '\.')[0,1] -join '.'
        $expectedMajorMinor = ($PYTHON_VERSION -split '\.')[0,1] -join '.'

        if ($majorMinor -eq $expectedMajorMinor) {
            Write-Success "Python version verified: $pyVersion"
            Write-Output "PYTHON_VERSION_OK=true"
            return $true
        } else {
            Write-Warning "Python version mismatch. Expected: $PYTHON_VERSION, Got: $pyVersion"
            Write-Info "Environment will need to be recreated..."
            Write-Output "PYTHON_VERSION_OK=false"
            return $false
        }
    } catch {
        Write-Error "Failed to verify Python version: $_"
        Write-Output "PYTHON_VERSION_OK=false"
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

    Write-Info "Step 3: Managing conda environment..."
    Write-Info "Using conda at: $CondaPath"

    if (-not (Test-CondaEnv $CondaPath)) {
        New-CondaEnv $CondaPath
    }

    if (-not (Test-PythonVersion $CondaPath)) {
        Write-Info "Recreating environment with correct Python version..."
        & "$CondaPath" remove -n $CONDA_ENV_NAME --all -y 2>$null
        New-CondaEnv $CondaPath
        Test-PythonVersion $CondaPath  # This should succeed now
    }
}

# Run main function with arguments
Main $args[0]