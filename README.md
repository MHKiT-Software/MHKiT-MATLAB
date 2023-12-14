# MHKiT-MATLAB

[![Linux and macOS Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml) [![Windows Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3928405.svg)](https://doi.org/10.5281/zenodo.3928405)

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

### Software Requirements

MHKiT-MATLAB utilizes Python functions from MHKiT-Python and requires the user to have
compatible versions of Python and MHKiT-Python installed.

MHKiT-MATLAB supports the following combinations of MATLAB and Python versions.[^1]

|      | R2021b | R2022a | R2022b | R2023a | R2023b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.11 |        |        |        |        | ✓      |
| 3.10 |        |        | ✓      | ✓      | ✓      |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8  | ✓      | ✓      | ✓      | ✓      |        |

Before installing MHKiT-MATLAB, please ensure your system has compatible versions of Python and MATLAB installed per the table above.

### Installation Guide

For complete installation instructions, please visit the [installation guide](https://mhkit-software.github.io/MHKiT/installation.html).

## Unit Tests

Our software undergoes rigorous testing through an [extensive suite of unit tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions), ensuring its stability and reliability. These tests simulate a user's machine, and failures on GitHub Actions may not necessarily indicate actual issues but could be artifacts of the build environment.

### Test Matrices

#### Linux (`ubuntu-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.11 | -      | -      | -      | -      | ✓      |
| 3.10 | -      | -      | ✓      | ✓      | ✓      |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8  | ✓      | ✓      | ✓      | ✓      | -      |

#### macOS (`macos-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.11 | -      | -      | -      | -      | X      |
| 3.10 | -      | -      | ✓      | ✓      | ✓      |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8  | ✓      | ✓      | ✓      | ✓      | ✓      |

#### Windows (`windows-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.11 | -      | -      | -      | -      | X      |
| 3.10 | -      | -      | X      | X      | X      |
| 3.9  | ✓      | ✓      | X      | X      | X      |
| 3.8  | ✓      | ✓      | X      | X      | X      |

### Legend

- ✓: GitHub Actions test passed.
- X: GitHub Actions test failed; consider using a tested version if issues arise.
- -: MATLAB/Python version not compatible.

If users face difficulties, using a tested version can help determine whether the problem is in the code or the GitHub Actions environment. Feel free to provide any additional details or specifications for further clarity.

## Copyright and License

MHKiT is copyright through the National Renewable Energy Laboratory,
Pacific Northwest National Laboratory, and Sandia National Laboratories.
The software is distributed under the Revised BSD License.

See [copyright and license](https://mhkit-software.github.io/MHKiT/license.html) for more information.

[^1]: For a comprehensive list of compatible MATLAB/Python versions, refer to the [MathWorks Python
    Compatibility Documentation](https://www.mathworks.com/support/requirements/python-compatibility.html).
