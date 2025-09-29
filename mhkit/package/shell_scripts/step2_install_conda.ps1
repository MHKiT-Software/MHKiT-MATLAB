# MHKiT-MATLAB Step 2: Install Miniconda
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
# Function to install miniconda
function Install-Conda {
    Write-Info "Installing miniconda..."

    $minicondaUrl = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
    $tempDir = [System.IO.Path]::GetTempPath()
    $installerPath = Join-Path $tempDir "miniconda.exe"

    # Download miniconda installer
    Write-Info "Downloading miniconda installer..."
    try {
        Invoke-WebRequest -Uri $minicondaUrl -OutFile $installerPath
    } catch {
        Write-Error "Failed to download miniconda installer: $_"
        Write-Output "CONDA_INSTALLED=false"
        exit 1
    }

    # Install miniconda silently
    Write-Info "Installing miniconda to $env:USERPROFILE\miniconda3..."
    try {
        $installArgs = @("/S", "/InstallationType=JustMe", "/AddToPath=0", "/RegisterPython=0", "/D=$env:USERPROFILE\miniconda3")
        Start-Process -FilePath $installerPath -ArgumentList $installArgs -Wait -NoNewWindow
    } catch {
        Write-Error "Failed to install miniconda: $_"
        Write-Output "CONDA_INSTALLED=false"
        exit 1
    }

    # Remove installer
    Remove-Item $installerPath -Force -ErrorAction SilentlyContinue

    Write-Success "Miniconda installed successfully"

    # Output conda path for next steps
    $condaPath = Join-Path $env:USERPROFILE "miniconda3\Scripts\conda.exe"
    Write-Output "CONDA_PATH=$condaPath"
}

function Main {
    Write-Info "Step 2: Installing miniconda..."
    Install-Conda
}

Main
