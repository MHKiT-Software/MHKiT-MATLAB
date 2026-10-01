<#
.SYNOPSIS
    Run the full MHKiT-MATLAB unit test suite, including tests that call
    MHKiT-Python, on Windows

.DESCRIPTION
    From the repository root, in PowerShell:

        scripts\run_python_tests_windows.ps1 C:\path\to\python.exe

    The Python interpreter must have mhkit and mhkit_python_utils installed:

        conda create -n mhkit_matlab -c conda-forge python=3.12 "numpy>=2" pip netcdf4 hdf5
        conda activate mhkit_matlab
        pip install "mhkit[all]==1.1.2" "pandas<3"  # pecos 1.0.0 check_delta fails with pandas 3
        pip install -e .
        scripts\run_python_tests_windows.ps1 (python -c "import sys; print(sys.executable)")

    If script execution is blocked, run once:
        Set-ExecutionPolicy -Scope CurrentUser RemoteSigned

    MATLAB is located from (in order): $env:MATLAB_EXE, `matlab` on PATH,
    the newest C:\Program Files\MATLAB\R20xx?\bin\matlab.exe.

    To isolate a users python config MATLAB is started with a minimal PATH
    and a throw-away preferences directory. This means any user `pyenv` MATLAB
    configuration is ignored.

    The Python environment folders are placed first on PATH so that the
    OutOfProcess Python host loads the environment's DLLs (expat, OpenSSL,
    HDF5, netCDF) instead of the versions that ship with MATLAB. Without this
    imports fail with errors like "DLL load failed while importing pyexpat".
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$PythonExe
)

$ErrorActionPreference = 'Stop'
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')
$PythonExe = (Resolve-Path $PythonExe).Path

function Find-Matlab {
    if ($env:MATLAB_EXE) { return $env:MATLAB_EXE }
    $onPath = Get-Command matlab.exe -ErrorAction SilentlyContinue
    if ($onPath) { return $onPath.Source }
    $installs = Get-ChildItem 'C:\Program Files\MATLAB' -Directory -Filter 'R20*' -ErrorAction SilentlyContinue |
        Sort-Object Name
    if ($installs) {
        return Join-Path $installs[-1].FullName 'bin\matlab.exe'
    }
    throw 'Could not find MATLAB. Set $env:MATLAB_EXE to the path of matlab.exe'
}

& $PythonExe -c "import mhkit, mhkit_python_utils; print('MHKiT-Python', mhkit.__version__)"
if ($LASTEXITCODE -ne 0) {
    throw "mhkit and mhkit_python_utils must be installed for $PythonExe"
}

# conda environments keep their DLLs in Library\bin, venvs next to python.exe
$envRoot = Split-Path $PythonExe -Parent
$envPath = @(
    $envRoot,
    (Join-Path $envRoot 'Library\mingw-w64\bin'),
    (Join-Path $envRoot 'Library\usr\bin'),
    (Join-Path $envRoot 'Library\bin'),
    (Join-Path $envRoot 'Scripts')
) | Where-Object { Test-Path $_ }

$matlabExe = Find-Matlab
$prefDir = Join-Path ([System.IO.Path]::GetTempPath()) ("mhkit_matlab_prefs_" + [System.Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $prefDir | Out-Null

Write-Host "MATLAB:   $matlabExe"
Write-Host "Python:   $PythonExe"
Write-Host "Prefdir:  $prefDir (temporary)"

$pythonForMatlab = $PythonExe.Replace("'", "''")
$matlabCommand = "addpath(fullfile('mhkit','tests')); results = runPythonTests('$pythonForMatlab'); disp(table(results)); assertSuccess(results);"

$savedPath = $env:PATH
$savedPrefDir = $env:MATLAB_PREFDIR
try {
    Push-Location $repoRoot
    $env:PATH = (($envPath + @("$env:SystemRoot\System32", "$env:SystemRoot")) -join ';')
    $env:MATLAB_PREFDIR = $prefDir
    # -wait keeps matlab.exe attached so the exit code is reported.
    & $matlabExe -batch $matlabCommand -wait
    $exitCode = $LASTEXITCODE
}
finally {
    Pop-Location
    $env:PATH = $savedPath
    if ($null -eq $savedPrefDir) { Remove-Item Env:\MATLAB_PREFDIR -ErrorAction SilentlyContinue } else { $env:MATLAB_PREFDIR = $savedPrefDir }
    Remove-Item -Recurse -Force $prefDir -ErrorAction SilentlyContinue
}
exit $exitCode
