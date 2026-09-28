# MHKiT-MATLAB Installation Guide

## Installation for Users

See: [MHKiT-MATLAB User Installation Guide](https://mhkit-software.github.io/MHKiT/matlab_installation.html)

## Installation for Developers

### Prerequisites

1. **Fork** the repository on GitHub: https://github.com/MHKiT-Software/MHKiT-MATLAB

2. **Clone** your fork to your local machine:

   ```bash
   git clone https://github.com/YOUR-USERNAME/MHKiT-MATLAB.git
   cd MHKiT-MATLAB
   ```

3. **Add upstream** as a remote:
   ```bash
   git remote add upstream https://github.com/MHKiT-Software/MHKiT-MATLAB.git
   ```

### Development Environment Setup

1. **Uninstall** the MHKiT toolbox if already installed:

   - Navigate to Home > Add-Ons > Manage Add-Ons
   - Right-click on "mhkit" > "Uninstall"

2. Create a conda environment for development:

```bash
conda create -n mhkit-matlab-dev python=311.
```

```bash
conda activate mhkit-matlab-dev
```

3. Install MHKiT-Python:

Note: Developers should be careful to install the correct version of MHKiT-Python. The context of
correct version depends on the specific task you are working on. For development against the release
version of MHKiT-MATLAB, install the latest stable release of MHKiT-Python from conda-forge:

```bash
conda install -c conda-forge mhkit
```

To install the latest development version of MHKiT-Python from source, clone the MHKiT-Python
repository:

```bash
git clone https://github.com/MHKiT-Software/MHKiT-Python.git
```

Optional, checkout the develop branch:

```bash
git checkout develop
```

Install MHKiT-Python in from source:

```bash
pip install -e .
```

Verify MHKiT-Python installation:

```bash
python -c "import mhkit; print(mhkit.__version__)"
```

4. Install mhkit-python-utils from source

Note: mhkit-python-utils is a required dependency for MHKiT-MATLAB. This python code handles any
MATLAB-Python interoperability functionality specific for MHKiT-MATLAB.

```bash
pip install -e .
```

5. Determine the path to the Python executable in your conda environment:

```bash
python -c "import sys; print(sys.executable)"
```

6. Add the mhkit conda environment to MATLAB:

In the MATLAB command window, run:

```matlab
pyenv(Version="<python executable location from above>")
```

Restart MATLAB

7. **Add to MHKiT-MATLAB mhkit and examples folder to the MATLAB path**:

```matlab
addpath(genpath('path/to/MHKiT-MATLAB/mhkit'))
```

If you encounter issues with Python integration, please [submit an issue](https://github.com/MHKiT-Software/MHKiT-MATLAB/issues/new).
