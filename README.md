# MHKiT-MATLAB

[![MATLAB Modules](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/matlab_only_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/matlab_only_unit_tests.yml) [![MATLAB Examples](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/matlab_only_examples.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/matlab_only_examples.yml) [![macOS Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml) [![Windows Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3928405.svg)](https://doi.org/10.5281/zenodo.3928405)

MHKiT-MATLAB is a MATLAB package designed for marine renewable energy applications to assist in
data processing and visualization. The software package include functionality for:

- Data processing
- Data visualization
- Data quality control
- Resource assessment
- Device performance
- Device loads

See the [documentation](https://mhkit-software.github.io/MHKiT/) for more information about MHKiT.

## Installation

### Quick Install

1. Download the MHKiT toolbox, `mhkit_v<version>.mltbx`, from the
   [latest release](https://github.com/MHKiT-Software/MHKiT-MATLAB/releases/latest).
2. Open the downloaded file in MATLAB (double-click it or drag it into the Command Window), or run:

   ```matlab
   matlab.addons.install("mhkit_v1.1.0.mltbx");
   ```

3. Verify the install:

   ```matlab
   matlab.addons.installedAddons
   ```

The acoustics, dolfyn, mooring, power, river, tidal, and most wave functions are native MATLAB and work after
this step. The loads, qc, river Delft3D, and some wave and utils functions call MHKiT-Python, which also requires
Python and MHKiT-Python, see [Software Requirements](#software-requirements) and the
[MHKiT MATLAB Installation Instructions](https://mhkit-software.github.io/MHKiT/matlab_installation.html).

To upgrade, install the new `.mltbx` over the existing one. To uninstall, go to
Home > Add-Ons > Manage Add-Ons, right-click on "Marine and Hydrokinetic Toolkit (MHKiT)", and select "Uninstall".

### Software Requirements

Some MHKiT-MATLAB modules utilize Python functions from MHKiT-Python and require the user to have
compatible versions of Python and MHKiT-Python installed.

MHKiT-MATLAB supports the following combinations of MATLAB and Python versions.[^1]

| Python | R2023b | R2024a | R2024b | R2025a | R2025b | R2026a | R2026b |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 3.13   | -      | -      | -      | -      | -      | ✓      | ✓      |
| 3.12   | -      | -      | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.11   | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.10   | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      |

- ✓: MATLAB/Python versions compatible
- `-`: MATLAB/Python versions not compatible

The minimum supported MATLAB release is R2023b. MHKiT-Python 1.1 requires Python 3.10 or newer. R2026b also
supports Python 3.14, which has not yet been tested with MHKiT-Python and MHKiT-MATLAB.

Before installing MHKiT-MATLAB, please ensure your system has compatible versions of Python and MATLAB installed per the table above.

### Installation Guide

For complete installation instructions, please visit the [installation guide](https://mhkit-software.github.io/MHKiT/installation.html).

## Unit Tests

To ensure software reliability and stability, MHKiT-MATLAB [runs a suite of unit tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions)
on GitHub Actions. These tests simulate a user's machine, but they are not perfect. Unit test failures on
GitHub Actions may not necessarily indicate actual issues but could be artifacts of the build environment.
Users should consider using a tested configuration if issues arise.

### Test Matrices

#### MATLAB-only modules (no Python)

Tests for the native MATLAB modules (acoustics, dolfyn, mooring, power, river, tidal, and most wave functions)
run without Python.

| OS                       | R2023b | R2024a | R2024b | R2025a | R2025b | R2026a | R2026b |
| ------------------------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| macOS (`macos-15`)       | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      |
| Windows (`windows-2025`) | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      | ✓      |

The MATLAB examples are additionally run on the oldest and newest supported releases, R2023b and R2026b, on
both operating systems.

#### Full test suite with MHKiT-Python

The complete test suite, including modules that call MHKiT-Python, runs on the latest release, R2026b, and on R2025b on each OS.

| OS                       | MATLAB | Python | MHKiT-Python |
| ------------------------ | ------ | ------ | ------------ |
| macOS (`macos-15`)       | R2025b | 3.12   | 1.1.2        |
| macOS (`macos-15`)       | R2026b | 3.12   | 1.1.2        |
| Windows (`windows-2025`) | R2025b | 3.12   | 1.1.2        |
| Windows (`windows-2025`) | R2026b | 3.12   | 1.1.2        |

Linux is not currently tested. Other MATLAB/Python combinations listed in
[Software Requirements](#software-requirements) are expected to work but are not exercised in CI.

### Legend

- ✓: Tested on GitHub Actions.

## Development Notes

### Contributions

We encourage contributions through pull requests. Please submit your contributions via pull requests on this repository.

### Local Development

#### Setup

1. Uninstall the MHKiT toolbox if already installed:

   - Navigate to Home > Add-Ons > Manage Add-Ons > right-click on "Marine and Hydrokinetic Toolkit (MHKiT)" > "Uninstall"

2. Clone or download the MHKiT-MATLAB source code. If contributing code, fork the repository and submit a pull request. GitHub provides details on the forking and pull request process [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests).

3. Install the latest Python versions of `mhkit` and `mhkit_python_utils`.

   - Navigate to the `MHKiT-MATLAB` directory:
     - Install `mhkit-python` with all module dependencies:
       - `pip install "mhkit[all]" "pandas<3"` (pecos, used by the qc module, does not yet support pandas 3)
     - Install `mhkit-python-utils`:
       - `pip install -e .`

4. Add the `MHKiT-MATLAB/mhkit` folder and its subfolders to your MATLAB path.

### Local Unit Testing

Ensure code integrity by running unit tests locally before pushing changes to GitHub.

To execute all unit tests, run `mhkit/tests/runTests.m`. Unit test results will display in the command window.

The same test runs used by GitHub Actions can be run from a terminal in the repository root:

- Native MATLAB modules only (no Python required):
  - macOS: `scripts/run_matlab_only_tests_macos.sh`
  - Windows: `scripts\run_matlab_only_tests_windows.ps1`
- Full test suite, including modules that call MHKiT-Python:
  - macOS: `scripts/run_python_tests_macos.sh /path/to/python`
  - Windows: `scripts\run_python_tests_windows.ps1 C:\path\to\python.exe`

### Code Coverage

Code coverage reports are automatically generated when running `mhkit/tests/runTests.m` (refer to [Local Unit Testing](#local-unit-testing)). The HTML report is written to `mhkit/tests/coverage_report`, which is not tracked in the repository.

## Copyright and License

MHKiT is copyright through the National Laboratory of the Rockies,
Pacific Northwest National Laboratory, and Sandia National Laboratories.
The software is distributed under the Revised BSD License.

See [copyright and license](https://mhkit-software.github.io/MHKiT/license.html) for more information.

[^1]:
    For a comprehensive list of compatible MATLAB/Python versions, refer to the [MathWorks Python
    Compatibility Documentation](https://www.mathworks.com/support/requirements/python-compatibility.html).
