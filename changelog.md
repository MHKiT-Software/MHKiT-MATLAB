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
