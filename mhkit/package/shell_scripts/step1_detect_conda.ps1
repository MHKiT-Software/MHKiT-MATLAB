# MHKiT-MATLAB Step 1: Conda Detection
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
# Function to check if conda is available
function Find-Conda {
    Write-Info "Starting comprehensive conda detection..."

    # Step 1: Check if conda is already in PATH and working
    try {
        $condaCmd = Get-Command conda -ErrorAction SilentlyContinue
        if ($condaCmd) {
            Write-Info "Found conda in PATH"
            if (Test-Conda $condaCmd.Source) {
                Write-Output "CONDA_DETECTED=true"
                Write-Output "CONDA_PATH=$($condaCmd.Source)"
                return
            } else {
                Write-Warning "Conda in PATH is not working properly"
            }
        }
    } catch {
        # Continue to other detection methods
    }

    # Step 2: Check environment variables
    $envValue = [Environment]::GetEnvironmentVariable("CONDA_EXE")
    if ($envValue) {
        Write-Info "Checking environment variable CONDA_EXE: $envValue"
        $condaFromEnv = ""
        if ("CONDA_EXE" -eq "CONDA_EXE") {
            $condaFromEnv = $envValue
        } elseif ("CONDA_EXE" -eq "CONDA_PREFIX") {
            $condaFromEnv = Join-Path $envValue "Scripts\conda.exe"
        }

        if ($condaFromEnv -and (Test-Conda $condaFromEnv)) {
            Write-Output "CONDA_DETECTED=true"
            Write-Output "CONDA_PATH=$condaFromEnv"
            return
        }
    }
    $envValue = [Environment]::GetEnvironmentVariable("CONDA_PREFIX")
    if ($envValue) {
        Write-Info "Checking environment variable CONDA_PREFIX: $envValue"
        $condaFromEnv = ""
        if ("CONDA_PREFIX" -eq "CONDA_EXE") {
            $condaFromEnv = $envValue
        } elseif ("CONDA_PREFIX" -eq "CONDA_PREFIX") {
            $condaFromEnv = Join-Path $envValue "Scripts\conda.exe"
        }

        if ($condaFromEnv -and (Test-Conda $condaFromEnv)) {
            Write-Output "CONDA_DETECTED=true"
            Write-Output "CONDA_PATH=$condaFromEnv"
            return
        }
    }

    # Step 3: Check all common installation paths systematically
    Write-Info "Checking common conda installation paths..."
    $pathsToCheck = @(
        [Environment]::ExpandEnvironmentVariables("%USERPROFILE%\miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("%USERPROFILE%\miniconda3\condabin\conda.bat"),
        [Environment]::ExpandEnvironmentVariables("%USERPROFILE%\Miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("%USERPROFILE%\Anaconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("%USERPROFILE%\anaconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("C:\ProgramData\Miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("C:\ProgramData\Anaconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("C:\tools\miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("C:\tools\Miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("%LOCALAPPDATA%\Continuum\miniconda3\Scripts\conda.exe"),
        [Environment]::ExpandEnvironmentVariables("%LOCALAPPDATA%\Continuum\anaconda3\Scripts\conda.exe")
    )

    foreach ($path in $pathsToCheck) {
        if (Test-Conda $path) {
            Write-Output "CONDA_DETECTED=true"
            Write-Output "CONDA_PATH=$path"
            return
        }
    }

    Write-Info "No working conda installation found"
    Write-Output "CONDA_DETECTED=false"
}

function Main {
    Write-Info "Step 1: Detecting conda installation..."
    Find-Conda
}

Main
