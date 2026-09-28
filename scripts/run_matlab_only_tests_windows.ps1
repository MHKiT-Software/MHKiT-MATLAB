<#
.SYNOPSIS
    Run the MATLAB-only test suites on Windows,

.DESCRIPTION
    From the repository root, in PowerShell:

        scripts\run_matlab_only_tests_windows.ps1            # unit tests (default)
        scripts\run_matlab_only_tests_windows.ps1 tests      # unit tests
        scripts\run_matlab_only_tests_windows.ps1 examples   # example live scripts
        scripts\run_matlab_only_tests_windows.ps1 all        # both

    If script execution is blocked, run once:
        Set-ExecutionPolicy -Scope CurrentUser RemoteSigned

    MATLAB is located from (in order): $env:MATLAB_EXE, `matlab` on PATH,
    the newest C:\Program Files\MATLAB\R20xx?\bin\matlab.exe.

    To isolate a users python config MATLAB is started with a minimal PATH
    and a throw-away preferences directory. This means any user `pyenv` MATLAB 
    configuration is ignored
#>
[CmdletBinding()]
param(
    [ValidateSet('tests', 'examples', 'all')]
    [string]$What = 'tests'
)

$ErrorActionPreference = 'Stop'
$repoRoot = Resolve-Path (Join-Path $PSScriptRoot '..')

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

switch ($What) {
    'tests'    { $cmd = 'results = runMatlabOnlyTests();' }
    'examples' { $cmd = 'results = runMatlabOnlyExamples();' }
    'all'      { $cmd = 'results = [runMatlabOnlyTests(), runMatlabOnlyExamples()];' }
}

$matlabExe = Find-Matlab
$prefDir = Join-Path ([System.IO.Path]::GetTempPath()) ("mhkit_matlab_prefs_" + [System.Guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $prefDir | Out-Null

Write-Host "MATLAB:   $matlabExe"
Write-Host "Running:  $What"
Write-Host "Prefdir:  $prefDir (temporary)"

$matlabCommand = "addpath(fullfile('mhkit','tests')); $cmd disp(table(results)); assertSuccess(results);"

$savedPath = $env:PATH
$savedPrefDir = $env:MATLAB_PREFDIR
try {
    Push-Location $repoRoot
    $env:PATH = "$env:SystemRoot\System32;$env:SystemRoot"
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
