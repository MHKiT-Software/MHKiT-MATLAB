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

See the [MHKiT MATLAB Installation Instructions](https://mhkit-software.github.io/MHKiT/matlab_installation.html).

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

To ensure software reliability and stability. MHKiT-MATLAB software [runs a suite of unit tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions) using the above MATLAB/Python compatibility matrix on Linux (`ubuntu-latest`), MacOS (`macos-latest`) and Windows (`windows-latest`). These tests simulate a user's machine, but they are not perfect. Unit test failures on GitHub Actions may not necessarily indicate actual issues but could be artifacts of the build environment. Users should consider using a tested version if issues arise.

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

|      | R2022b | R2023a | R2023b | R2024a | R2024b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.12 | -      | -      | -      | -      | ✓      |
| 3.11 | -      | -      | ✓      | ✓      | ✓      |
| 3.10 | ✓      | ✓      | ✓      | ✓      | ✓      |

### Legend

- ✓: GitHub Actions unit test passed.
- X: GitHub Actions unit test failed; consider using a tested version if issues arise.
- `-`: MATLAB/Python version not compatible.

## Development Notes

### Contributions

We encourage contributions through pull requests. Please submit your contributions via pull requests on this repository.

### Local Development

#### Setup

1. Uninstall the MHKiT toolbox if already installed:

   - Navigate to Home > Add-Ons > Manage Add-Ons > right-click on "mhkit" > "Uninstall"

2. Clone or download the MHKiT-MATLAB source code. If contributing code, fork the repository and submit a pull request. GitHub provides details on the forking and pull request process [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests).

3. Install the latest Python versions of `mhkit` and `mhkit_python_utils`.

   - Navigate to the `MHKiT-MATLAB` directory:
     - Install `mhkit-python`:
       - `pip install mhkit`
     - Install `mhkit-python-utils`:
       - `pip install -e .`

4. Add the `MHKiT-MATLAB/mhkit` folder and its subfolders to your MATLAB path.

### Local Unit Testing

Ensure code integrity by running unit tests locally before pushing changes to GitHub.

To execute all unit tests, run `mhkit/tests/runTests.m`. Unit test results will display in the command window.

### Code Coverage

Code coverage reports are automatically generated during unit testing (refer to [Local Unit Testing](#local-unit-testing) for execution instructions). `mhkit/tests/coverage_report` contains results from the most recent code coverage report.

## Copyright and License

MHKiT is copyright through the National Laboratory of the Rockies,
Pacific Northwest National Laboratory, and Sandia National Laboratories.
The software is distributed under the Revised BSD License.

See [copyright and license](https://mhkit-software.github.io/MHKiT/license.html) for more information.

[^1]:
    For a comprehensive list of compatible MATLAB/Python versions, refer to the [MathWorks Python
    Compatibility Documentation](https://www.mathworks.com/support/requirements/python-compatibility.html).
