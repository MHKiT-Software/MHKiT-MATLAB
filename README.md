MHKiT-MATLAB
===================================
[![MHKiT-MATLAB Unix Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/unix_unit_tests.yml) [![MHKiT-MATLAB Windows Unit Tests](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml/badge.svg)](https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/windows_unit_tests.yml) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3928405.svg)](https://doi.org/10.5281/zenodo.3928405)

MHKiT-MATLAB is a MATLAB package designed for marine renewable energy applications to assist in
data processing and visualization.  The software package include functionality for:

* Quality control analysis
* Wave resource and wave performance
* River resource and device
* Tidal resource and device
* Data visualization
* Mechanical loads analysis
* Power quality assessments

See the [documentation](https://mhkit-software.github.io/MHKiT/) for more information about MHKiT.

Installation
------------------------

MHKiT-MATLAB consists of functions that wrap the MHKiT-Python functions. Therefore, MHKiT-MATLAB has specific Python package dependencies.

Before proceeding with the installation, please ensure compatibility with your MATLAB version and the recommended Python version combinations. Refer to the [Compatible Python and MATLAB Versions](#compatible-python-and-matlab-versions) section for the compatibility matrix.

For detailed installation instructions, please visit the [installation guide](https://mhkit-software.github.io/MHKiT/installation.html).

Compatible Python and MATLAB Versions
------------------------

MHKiT-MATLAB supports the following combinations of MATLAB and Python versions. For a comprehensive list of compatible MATLAB/Python versions, refer to [MathWorks Python Compatibility](https://www.mathworks.com/support/requirements/python-compatibility.html).

|       | R2021b | R2022a | R2022b | R2023a | R2023b |
| ----- | ------ | ------ | ------ | ------ | ------ |
| 3.11  |        |        |        |        | ✓      |
| 3.10  |        |        | ✓      | ✓      | ✓      |
| 3.9   | ✓      | ✓      | ✓      | ✓      | ✓      |
| 3.8   | ✓      | ✓      | ✓      | ✓      |        |


Copyright and license
------------------------
MHKiT is copyright through the National Renewable Energy Laboratory,
Pacific Northwest National Laboratory, and Sandia National Laboratories.
The software is distributed under the Revised BSD License.

See [copyright and license](https://mhkit-software.github.io/MHKiT/license.html) for more information.
