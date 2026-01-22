# MHKiT-MATLAB Step 4: Install Dependencies
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
# Function to install pre-install dependencies
function Install-PreDependencies {
    param([string]$CondaPath)

    Write-Info "Installing pre-install dependencies..."

    Write-Info "Executing: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf netcdf4 pandas numpy -c conda-forge -y"
    try {
        conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf netcdf4 pandas numpy -c conda-forge -y
    } catch {
        Write-Error "Failed to execute: conda run -n mhkit-matlab-env conda install pip hdf5 libnetcdf netcdf4 pandas numpy -c conda-forge -y"
        Write-Output "PRE_DEPS_INSTALLED=false"
        exit 1
    }

    Write-Success "Pre-install dependencies installed"
    Write-Output "PRE_DEPS_INSTALLED=true"
}

# Function to check if mhkit is already installed with correct version
function Test-MhkitInstalled {
    param([string]$CondaPath)

    Write-Info "Checking if mhkit-python is already installed..."

    try {
        $installedVersion = & "$CondaPath" run -n $CONDA_ENV_NAME python -c "import mhkit; print(mhkit.__version__)" 2>$null
        if ($installedVersion -like "*$MHKIT_VERSION*") {
            Write-Success "mhkit-python v$MHKIT_VERSION already installed"
            Write-Output "MHKIT_INSTALLED=true"
            return $true
        } else {
            Write-Info "mhkit-python not installed or wrong version (found: $installedVersion)"
            Write-Output "MHKIT_INSTALLED=false"
            return $false
        }
    } catch {
        Write-Info "mhkit-python not installed or wrong version"
        Write-Output "MHKIT_INSTALLED=false"
        return $false
    }
}

# Function to install mhkit-python
function Install-Mhkit {
    param([string]$CondaPath)

    Write-Info "Installing mhkit-python v$MHKIT_VERSION..."

    try {
        & "$CondaPath" run -n $CONDA_ENV_NAME conda install -c conda-forge mhkit=$MHKIT_VERSION -y

        Write-Success "mhkit-python installed successfully"
        Write-Output "MHKIT_INSTALLED=true"

        # Verify installation
        Write-Info "Verifying mhkit-python installation..."
        $installedVersion = & "$CondaPath" run -n $CONDA_ENV_NAME python -c "import mhkit; print(mhkit.__version__)"

        if ($installedVersion -like "*$MHKIT_VERSION*") {
            Write-Success "Installation verified: $installedVersion"
            Write-Output "MHKIT_VERIFIED=true"
        } else {
            Write-Error "Version verification failed. Expected: v$MHKIT_VERSION, Got: $installedVersion"
            Write-Output "MHKIT_VERIFIED=false"
            exit 1
        }
    } catch {
        Write-Error "Failed to install mhkit-python: $_"
        Write-Output "MHKIT_INSTALLED=false"
        exit 1
    }
}

# Function to install mhkit_python_utils
function Install-MhkitPythonUtils {
    param([string]$CondaPath)

    Write-Info "Installing mhkit_python_utils..."

    # Create temporary directory for download
    $tempDir = [System.IO.Path]::GetTempPath()
    $utilsTempDir = Join-Path $tempDir "mhkit_utils_$(Get-Random)"
    New-Item -ItemType Directory -Path $utilsTempDir -Force | Out-Null

    $utilsUrl = "https://github.com/MHKiT-Software/MHKiT-MATLAB/releases/download/v0.6.0/mhkit_python_utils_v0.6.0.zip"
    $zipPath = Join-Path $utilsTempDir "mhkit_python_utils.zip"

    Write-Info "Downloading mhkit_python_utils from: $utilsUrl"

    try {
        # Download the zip file
        Invoke-WebRequest -Uri $utilsUrl -OutFile $zipPath

        # Extract the package
        Write-Info "Extracting mhkit_python_utils..."
        Expand-Archive -Path $zipPath -DestinationPath $utilsTempDir -Force

        # Find the extracted directory
        $utilsDir = Get-ChildItem -Path $utilsTempDir -Directory | Where-Object { $_.Name -like "*mhkit_python_utils*" } | Select-Object -First 1

        if (-not $utilsDir) {
            Write-Error "Could not find extracted mhkit_python_utils directory"
            Write-Output "UTILS_INSTALLED=false"
            exit 1
        }

        # Install in development mode
        Write-Info "Installing mhkit_python_utils with pip install -e ..."
        Set-Location $utilsDir.FullName
        & "$CondaPath" run -n $CONDA_ENV_NAME pip install -e .

        Write-Success "mhkit_python_utils installed successfully"
        Write-Output "UTILS_INSTALLED=true"

        # Verify installation
        Write-Info "Verifying mhkit_python_utils installation..."
        try {
            $utilsVersion = & "$CondaPath" run -n $CONDA_ENV_NAME python -c "import mhkit_python_utils; print(mhkit_python_utils.__version__)" 2>$null
            Write-Info "mhkit_python_utils version: $utilsVersion"
        } catch {
            Write-Info "mhkit_python_utils version: unknown"
        }

        # Clean up temporary directory
        Remove-Item $utilsTempDir -Recurse -Force -ErrorAction SilentlyContinue

    } catch {
        Write-Error "Failed to install mhkit_python_utils: $_"
        Write-Output "UTILS_INSTALLED=false"
        Remove-Item $utilsTempDir -Recurse -Force -ErrorAction SilentlyContinue
        exit 1
    }
}

function Main {
    param([string]$CondaPath)

    if (-not $CondaPath) {
        Write-Error "No conda path provided. Usage: script.ps1 <conda_path>"
        exit 1
    }

    Write-Info "Step 4: Installing dependencies..."
    Write-Info "Using conda at: $CondaPath"

    Install-PreDependencies $CondaPath

    if (-not (Test-MhkitInstalled $CondaPath)) {
        Install-Mhkit $CondaPath
    }

    Install-MhkitPythonUtils $CondaPath
}

# Run main function with arguments
Main $args[0]
