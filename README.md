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
| 3.11 | -      | -      | -      | -      | ✓      |
| 3.10 | -      | -      | ✓      | ✓      | ✓      |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8  | ✓      | ✓      | ✓      | ✓      | -      |

- ✓: MATLAB/Python versions compatible
- `-`: MATLAB/Python versions not compatible

Before installing MHKiT-MATLAB, please ensure your system has compatible versions of Python and MATLAB installed per the table above.

### Installation Guide

For complete installation instructions, please visit the [installation guide](https://mhkit-software.github.io/MHKiT/installation.html).

## Unit Tests

To ensure software reliability and stability. MHKiT-MATLAB software [runs a suite of unit tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions) using the above MATLAB/Python compatibility matrix on Linux (`ubuntu-latest`), MacOS (`macos-latest`) and Windows (`windows-latest`). These tests simulate a user's machine, but they are not perfect. Unit test failures on GitHub Actions may not necessarily indicate actual issues but could be artifacts of the build environment. Users should consider using a tested version if issues arise.

### Test Matrices

The test matrices below detail the current state of unit testing. An "X" indicates a OS/MATLAB/Python version with a failing MHKiT-MATLAB unit test on GitHub Action that is due to Actions environment.

#### Linux (`ubuntu-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b |
| ---- | ------ | ------ | ------ | ------ | ------ |
| 3.11 | -      | -      | -      | -      | ✓      |
| 3.10 | -      | -      | ✓      | ✓      | ✓      |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8  | ✓      | ✓      | ✓      | ✓      | -      |

#### macOS (`macos-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b                     |
| ---- | ------ | ------ | ------ | ------ | -------------------------- |
| 3.11 | -      | -      | -      | -      | X[^macos-311-action-error] |
| 3.10 | -      | -      | ✓      | ✓      | ✓                          |
| 3.9  | ✓      | ✓      | ✓      | ✓      | ✓                          |
| 3.8  | ✓      | ✓      | ✓      | ✓      | ✓                          |

#### Windows (`windows-latest`)

|      | R2021b | R2022a | R2022b | R2023a | R2023b                            |
| ---- | ------ | ------ | ------ | ------ | --------------------------------- |
| 3.11 | -      | -      | -      | -      | X[^windows-gt-2022a-action-error] |
| 3.10 | -      | -      | X      | X      | X                                 |
| 3.9  | ✓      | ✓      | X      | X      | X                                 |
| 3.8  | ✓      | ✓      | X      | X      | X                                 |

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

MHKiT is copyright through the National Renewable Energy Laboratory,
Pacific Northwest National Laboratory, and Sandia National Laboratories.
The software is distributed under the Revised BSD License.

See [copyright and license](https://mhkit-software.github.io/MHKiT/license.html) for more information.

[^1]:
    For a comprehensive list of compatible MATLAB/Python versions, refer to the [MathWorks Python
    Compatibility Documentation](https://www.mathworks.com/support/requirements/python-compatibility.html).

[^macos-311-action-error]:

[macOS Python 3.11 action error](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/runs/7143629923/job/19455784365#step:20:3373)

    ```
    Python Error: ImportError:
          dlopen(/usr/local/miniconda/envs/MHKIT_CONDA_ENV/lib/python3.11/site-packages/netCDF4/_netCDF4.cpython-311-darwin.so,
          0x0002): Symbol not found: (_SSL_get0_group_name)
            Referenced from:
            '/usr/local/miniconda/envs/MHKIT_CONDA_ENV/lib/libcurl.4.dylib'
            Expected in:
            '/Users/runner/hostedtoolcache/MATLAB/2023.2.999/x64/MATLAB.app/bin/maci64/libssl.3.dylib'
    ```

[^windows-gt-2022a-action-error]:
    [Windows MATLAB > 2022a action error](https://github.com/MHKiT-Software/MHKiT-MATLAB/pull/108/commits/63fbf91543392e2fc673a67af52e152ff15f2980)

    ```
    Error: Caught "std::exception" Exception message is:
      Unable to launch Simple server: Unable to launch C:\Program Files\MATLAB\R2023b\interprocess\bin\win64\pycli\MATLABPyHost.exe
      because: Peer process exited before transport type handshake. Exit status: 3221226505

      ERROR: MATLAB error Exit Status: 0x00000001
      exit status 1
    ```
