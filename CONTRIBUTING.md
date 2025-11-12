# Contributing to MHKiT-MATLAB

Thank you for your interest in contributing to MHKiT-MATLAB! This guide will help you get started with contributing code, documentation, and bug fixes.

## Quickstart

Here's the high-level contribution workflow:

1. **Fork** the MHKiT-MATLAB repository
2. **Set up** your development environment (Python + MATLAB)
3. **Create** a feature branch from `develop`
4. **Write** your code with proper docstrings and tests
5. **Test locally** using the docstring checker and unit tests
6. **Add** an entry to `changelog.md`
7. **Push** your changes and open a pull request to `develop`
8. **CI runs** automated tests across all supported MATLAB/Python versions
9. **Address** review feedback from maintainers
10. **Merge** into `develop` once approved

## Fork and Setup

### Fork the Repository

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

3. **Install Python dependencies**:

Install mhkit-python-utils from source

```bash
pip install -e .
```

4. Determine the path to the Python executable in your conda environment:

```bash
python -c "import sys; print(sys.executable)"
```

4. Add the mhkit conda environment to MATLAB:

In the MATLAB command window, run:

```matlab
pyenv(Version="<python executable location from above>")
```

Restart MATLAB

5. **Add to MHKiT-MATLAB mhkit and examples folder to the MATLAB path**:

```matlab
addpath(genpath('path/to/MHKiT-MATLAB/mhkit'))
```

If you encounter issues with Python integration, please sumbit an issue.

## Development Workflow

### Branching Strategy

- **`master`** - Production releases only
- **`develop`** - Main development branch ← **target your PRs here**

### Creating a Feature Branch

1. **Sync with upstream**:

   ```bash
   git checkout develop
   git fetch upstream
   git merge upstream/develop
   ```

2. **Create your feature branch**:

   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make changes and commit regularly**:

   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

4. **Keep your branch updated**:
   ```bash
   git fetch upstream
   git rebase upstream/develop
   ```

## Writing Code

## Documentation

### Docstring Format

All public functions **must** have properly formatted docstrings. This is loosely enforced by automated CI checks. Where a specific format is not specified, developers should use the template below as a guide.

### Docstring Function Template

```matlab
function result = descriptive_function_name(input_1, input_2, input_3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Brief one-line description of what the function does
%
% Optional extended description paragraph. Can span multiple lines
% to provide additional context, mathematical background, data sources,
% or important behavioral notes about the function.
%
% Parameters
% ------------
% input_1 : type [units]
%   Description of first parameter. Can span multiple lines.
%   Additional details aligned with first line. Use two spaces for indentation
%     input_1.fieldname_1 : type [units]
%       Description of field
%     input_1.fieldname_2 : type [units]
%       Description of field
% input_2 : type [units]
%   Description of second parameter
% input_t : type [units] (optional)
%   Description of second parameter. Mark optional parameters
%   with (optional) after the type.
%
% Returns
% ---------
% output : type
%   Description of return value. For structure outputs,
%   document nested fields with indentation:
%     output.field1 : type [units]
%       Field 1 description
%     output.field2 : type [units]
%       Field description
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    input_1 struct
    input_2 {mustBeInteger, mustBePositive}
    input_3 {mustBeInteger} = 2
end

arguments (Output)
    result struct
end

% Function implementation goes here

end
```

### Use the `arguments` Block

All functions must use [MATLAB's `arguments` block](https://www.mathworks.com/help/matlab/ref/arguments.html) for input/output type validation.

#### Simple Example

```matlab
function result = calculate_energy(power, duration)

arguments
    power double {mustBeNumeric, mustBeFinite}
    duration double {mustBePositive}
end

% Function implementation
result = power * duration;

end
```

#### Complex Example (Struct Inputs/Outputs)

```matlab
function mspl = band_sound_pressure_level(spsd, octave, base, fmin, fmax)

arguments (Input)
    spsd struct
    octave {mustBeInteger, mustBePositive}
    base {mustBeInteger} = 2
    fmin {mustBeNumeric, mustBePositive} = 10
    fmax {mustBeNumeric, mustBePositive} = 100000
end

arguments (Output)
    mspl struct
end

% Function implementation
mspl = struct();

end
```

#### Common Validation Functions

- `mustBeNumeric` - Must be numeric type
- `mustBeFinite` - No NaN or Inf values
- `mustBePositive` - Greater than zero
- `mustBeInteger` - Integer values only
- `mustBeNonzero` - Cannot be zero
- `mustBeNonempty` - Cannot be empty

For custom validation, create functions in `mhkit/utils/` or module-specific directories.

### Checking Your Docstrings

Before submitting your PR, there is a python script to validate docstrings locally:

```bash
# Check all mhkit files
python scripts/check_docstrings.py

# Check specific directory
python scripts/check_docstrings.py --path mhkit/wave

# Verbose output
python scripts/check_docstrings.py --verbose
```

**Note**: The docstring checker actions run automatically on all PRs. Errors must be fixed before
merging to develop and master.

## Testing

### Running Tests Locally

Before submitting your PR, run the full test suite:

```matlab
% Navigate to tests directory
cd mhkit/tests

% Run all tests
runTests
```

### Writing Unit Tests

When adding new functionality, create tests in `mhkit/tests/`:

Use the patterns in the existing test files as a guide.

**Test coverage should include:**

- Normal operation
- Edge cases
- Error conditions
- Optional parameters

## Changelog

All contributions **must** include a changelog entry.

### Adding a Changelog Entry

1. Open `changelog.md`

2. Add your entry at the **top** under a new "Unreleased" section (if it doesn't exist):

```markdown
# Unreleased

## PR #XXX - Brief Feature Description

- Author: @your-github-username
- Description of changes
- Key additions:
  - Bullet point 1
  - Bullet point 2

[rest of existing changelog...]
```

### Changelog Format Guidelines

- Use format: `## PR #XXX - Title`
- Include your GitHub username
- List key changes as bullet points
- Reference related issues if applicable
- Place at the top of the file

**Example:**

```markdown
# Unreleased

## PR #185 - Add Wave Spectral Analysis Functions

- Author: @contributor
- Added native MATLAB implementations for wave spectral analysis
- Key additions:
  - `wave/resource/wave_spectral_moment.m` - Calculate spectral moments
  - `wave/resource/wave_spectral_bandwidth.m` - Calculate spectral bandwidth
  - Unit tests and example live script
- Closes #180
```

---

## Submitting Your PR

### Pull Request Process

1. **Push your changes** to your fork:

   ```bash
   git push origin feature/your-feature-name
   ```

2. **Create a pull request**:

   - Go to https://github.com/MHKiT-Software/MHKiT-MATLAB
   - Click "New Pull Request"
   - **Base branch**: `develop` ← Important!
   - **Compare branch**: your feature branch
   - Provide a clear description of your changes

3. **PR requirements**:

   - Target `develop` branch (NOT `master`)
   - All CI checks must pass:
     - Docstring validation
     - Unit tests (Linux, macOS, Windows)
     - Code compatibility checks
   - Changelog entry added
   - Tests added/updated for new functionality

4. **After submission**:
   - CI will automatically test across all supported MATLAB/Python versions
   - Address any review feedback
   - Once approved, maintainers will merge to `develop`

Note: Sometimes github actions fail for non code related reasons. If this happens, please reach
out in the discussions or issues page. It is a reasonable procedure to re-run failed actions, and
possibly disable failing combinations of matlab/python versions if they are not critical.

### PR Title Format

Use conventional commit format:

- `Feature: Add wave spectral bandwidth calculation`
- `<module>: Correct NDBC data parsing for missing values`
- `Actions: Fix action related bug`
- `Examples: Fix typo in acoustics example`
- `Fix: Incorrect sign in energy calculation`

## Code Standards

- **Style**: Use the templates provided above to ensure functions follow consistent style
- **Clarity**: Use meaningful variable names with units, i.e. `significant_wave_height_meters` over shorter variables like `hm0`.
- **Comments**: Use comments to explain complex logic, but avoid obvious comments
- **Testing**: Include unit tests

## Getting Help

- **Documentation**: https://mhkit-software.github.io/MHKiT/
- **Issues**: https://github.com/MHKiT-Software/MHKiT-MATLAB/issues
- **Discussions**: https://github.com/MHKiT-Software/MHKiT-MATLAB/discussions

## License

By contributing, you agree that your contributions will be licensed under the Revised BSD License. See [LICENSE](https://mhkit-software.github.io/MHKiT/license.html) for details.
