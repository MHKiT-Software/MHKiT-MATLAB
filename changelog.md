# Version 1.0.1

## Release Highlights

Bug fixes and documentation improvements

## Details

### Bug Fixes

- Issue #171 - Failing Loads: bin_statistics test
  - Fixed by @hivanov-nrel in PR #188
- MHKiT-Documentation build errors:
  <https://github.com/MHKiT-Software/MHKiT/pull/92#issuecomment-3457331254> - Fixed by @simmsa in PR #187

### Documentation Improvements

- Added GitHub Action to check docstring formatting
  - Verifies docstrings use a format that is compatible with MHKiT-Documentation build system
  - <https://github.com/MHKiT-Software/MHKiT-MATLAB/actions/workflows/check_docstrings.yml>
  - Implemented by @simmsa in PR #187
- Added CONTRIBUTING.md
  - Details guidelines for developers contributing to MHKiT-MATLAB
  - Implemented by @simmsa in PR #187
- Added INSTALL.md
  - Provides developer installation instructions for MHKiT-MATLAB
  - Implemented by @simmsa in PR #187
- Updated National Laboratory of the Rockies (NLR) branding
  - Implemented by @simmsa in PR #189

# Version 1.0.0

## Release Highlights

- New acoustics, and mooring modules
- Improvements to Wave, WDRT, and DOLFYN modules
- New examples for acoustics, mooring, and WEC-sim, and improved examples for DOLFYN and wave modules
- Multiple bug fixes and performance improvements

## PR #174 - Acoustics Module

- Authors: @hivanov-nrel, @simmsa
- Addition of acoustics module:
  - Reading and standardization of of OceanSonics icListen and OceanInstruments Soundtrap hydrophone files
  - Implementation of numerical computation sections IEC 62600-40 "Acoustic characterization of marine energy converters" standard
  - Spectrogram data visualizations using `plot_spectrogram`
  - Sound exposure level calculations with auditory weighting functions for 5 groups of marine mammals
  - Complete example livescript in examples/acoustics_example.mlx

## PR 176 - Wave Module Native MATLAB Implementation

- Authors: @MShabara, @simmsa
- Convert wave module functions to native MATLAB code:
  - wave/resource/standardize_wave_spectra_frequency.m:
  - wave/resource/frequency_moment.m:
  - wave/resource/energy_period.m:
  - wave/resource/significant_wave_height.m:
  - wave/resource/jonswap_spectrum.m:
  - wave/resource/pierson_moskowitz_spectrum.m:
  - wave/resource/surface_elevation.m:

## PR 153 - Mooring Module

- Author: @hivanov-nrel
- Addition of the mooring module, which includes
  - Read in MoorDyn files into MATLAB
  - Function to calculate mooring line lay length
  - Functions to visualize mooring line dynamics in 2D or 3D.
  - Example LiveScript demonstrating the functionality of the mooring module

## PR #173 - Dolfyn Turbulence Functionality

- Author: @simmsa, @browniea
- Addition of turbulence calculations for acoustic doppler instruments
  - Turbulence intensity, noise, and reynolds stress calculations
  - Updates to examples/adcp_example.mlx live script

## PR 170 - WDRT leftovers

- Author: @hivanov-nrel
- Expanded WDRT functionality by adding `automatic_hs_threshold` and `return_year_value` functions

## PR 169 - Code Compatibility Improvements

- Author: @simmsa
- Fix code compatibility issues detailed in issues #115, #116, #117, #118, #119, #120, #121, and #122

## Bug Fixes

- Issue #172 - Fix Failing MLER Test - @hivanov-nrel
- Issue #152 - Fix dimensionality differences in environmental_contours_example - @simmsa
- Issue #146, #74 - Finish WDRT Parity - @hivanov-nrel
- Issue #145 - Fix build errors in documentation - @simmsa
- Issue #114, #115, #116, #117, #118, #119, #121, #122 - MATLAB code compatibility issues

# MHKiT-MATLAB v0.6.0

## Release Highlights

- Added module and example for upcrossing analysis
- Expanded WDRT and extreme wave analysis capability
- Expanded DOLFYN functionality with performance optimizations and enhanced visualization
- Enhanced MATLAB-Python interoperability with robust type conversion
- Added support for MATLAB 2024b, Python 3.10-3.12, and MHKiT-Python v0.9.0
- Improved version control capabilities with diffable example formats

## Breaking Changes

- Removed support for:
  - MHKiT-Python versions below v0.9.0
  - MATLAB versions 2021b and 2022b
  - Python versions 3.8 and 3.9

Impact:

- Users must upgrade to MATLAB 2023a or newer
- Python environment must use version 3.10-3.12
- Users must upgrade to MHKiT-Python v0.9.0 or newer

## New Features

### Wave Module

- Added Upcrossing Analysis Functions and Example (@MShabara) #151
  - Enables time-series analysis of wave heights and periods using upcrossing methods
  - Documentation: [Upcrossing Analysis](https://github.com/MHKiT-Software/MHKiT-MATLAB/tree/master/mhkit/utils/upcrossing)
- Enhanced Environmental Contours (@hivanov-nrel) #147
  - Added support for IFORM and direct sampling methods
  - Improved extreme condition analysis capabilities
  - See examples in `environmental_contours_example.mlx`

### DOLFYN Module

- Enhanced Data Processing (@simmsa) #141
  - Added statistical functions:
    - Dimensional averaging [5b6d20a]
    - Speed and velocity calculations [8ad973b]
    - Time series analysis capabilities
  - Performance Optimizations:
    - Added caching for improved read_signature speed [85bf6e0]
    - Enhanced NetCDF attribute handling
  - Expanded visualization capabilities [8e4d81d]:
    - New subplot functionality
    - Custom colormaps similar to matplotlib (viridis, bluewhitered)
    - Histogram generation
    - Improved attribute visualization and customization
  - New comprehensive ADCP example in `adcp_example.mlx`

## Technical Improvements

### Compatibility

- Added official support for MATLAB 2024b
- Added support for MHKiT-Python v0.9.0 (@simmsa) #149
  - Compatible with Python 3.10-3.12
  - Enhanced type conversion system between MATLAB and Python

### MATLAB-Python Interoperability

- Added robust type conversion system:
  - New `typecast_spectra_to_mhkit_python`: Standardizes spectra data sent to Python
  - New `typecast_from_mhkit_python`: Standardizes Python output for MATLAB
    - Returns consistent struct format:
      - `data`: Native MATLAB type (scalar, array, or struct)
      - `type`: Data classification
      - `index`: Index information
    - Supports:
      - Scalar values
      - Python/NumPy arrays
      - pandas DataFrame/Series

### Version Control Improvements

- Added Diffable Examples (@simmsa) #157
  - MATLAB Markup (.m) versions of examples for better version tracking
  - Maintains Livescript functionality while enabling git-based collaboration
  - Improves documentation clarity and maintainability

## New Examples

- Strain Analysis Example (@simmsa, @akeeste) #156
  - Demonstrates structural load analysis techniques
  - Shows data processing workflow for strain measurements
  - Includes practical visualizations and analysis methods
- Upcrossing Analysis Example (@MShabara) #151
  - Shows wave height and period analysis workflow
  - Includes practical applications of statistical methods
- ADCP Data Processing Example (@simmsa) #141
  - Illustrates acoustic doppler current profiler data analysis
  - Demonstrates new DOLFYN visualization capabilities
  - Includes comprehensive processing workflow

## Bug Fixes

- Fixed WPTO hindcast download issue specific to `omni-directional_wave_power` [8284aa2]
  - [Issue #143](https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/143)
  - [Fixed in PR #144](https://github.com/MHKiT-Software/MHKiT-MATLAB/pull/144)
  - Improved robustness and clarity of hindcast downloads
- Corrected time scale representation in short-term extremes example (@MShabara) [#150, 3acf336]
- Fixed Delft3D masked array type conversion issues [2a290e5]
- Enhanced NetCDF attribute handling and extraction [f63002d]

## Dependencies

- MATLAB ≥ 2023a
- MHKiT-Python ≥ v0.9.0
- Python 3.10-3.12

## Contributors

Special thanks to all contributors who made this release possible:

- @hivanov-nrel
- @MShabara
- @simmsa
- @akeeste
- @rpauly18

# v0.5.0

## New Features

- Addition of the WEC Design Response Toolbox within the `wave` module #127
  - Estimate extreme sea states based on short term data
- Addition of the Delft3D input and analysis within the `tidal` module #124
  - Analyze modeled river/tidal flow data using same tools as ADCP and resource data

## Improvements

- More detailed and complete [installation instructions](https://mhkit-software.github.io/MHKiT/matlab_installation.html)
- Update MATLAB/Python compatibility matrix

## Fixes

- Allow user to specify surface elevation generation method #126
- Properly map the gamma parameter in the `jonswap` function #136
