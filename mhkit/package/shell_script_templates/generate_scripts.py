#!/usr/bin/env python3
"""
Script to generate shell installation scripts from Jinja2 templates and spec.json

This script reads the MHKiT specification from spec.json and generates
platform-specific shell scripts for Python environment setup.

Usage:
    python generate_scripts.py

Generated files:
    - ../scripts/install_mhkit_python_unix.sh (for macOS/Linux)
    - ../scripts/install_mhkit_python_windows.ps1 (for Windows)
"""

import json
import os
import sys
from pathlib import Path

try:
    from jinja2 import Environment, FileSystemLoader
except ImportError:
    print("Error: jinja2 is required to generate scripts")
    print("Install with: pip install jinja2")
    sys.exit(1)


def load_spec():
    """Load the MHKiT specification from spec.json"""
    script_dir = Path(__file__).parent
    spec_path = script_dir.parent / "+mhkit" / "spec.json"

    if not spec_path.exists():
        raise FileNotFoundError(f"spec.json not found at {spec_path}")

    with open(spec_path, 'r') as f:
        return json.load(f)


def setup_jinja_environment():
    """Set up Jinja2 environment with templates directory"""
    script_dir = Path(__file__).parent
    template_dir = script_dir

    return Environment(
        loader=FileSystemLoader(template_dir),
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True
    )


def replace_placeholders_in_string(text, replacements):
    """Replace all placeholders in a string with their values"""
    if not isinstance(text, str):
        return text

    result = text
    for placeholder, value in replacements.items():
        result = result.replace(f"<{placeholder}>", str(value))

    return result


def replace_placeholders_recursively(obj, replacements):
    """Recursively replace placeholders in nested structures"""
    if isinstance(obj, str):
        return replace_placeholders_in_string(obj, replacements)
    elif isinstance(obj, list):
        return [replace_placeholders_recursively(item, replacements) for item in obj]
    elif isinstance(obj, dict):
        return {key: replace_placeholders_recursively(value, replacements) for key, value in obj.items()}
    else:
        return obj


def preprocess_spec_for_templates(spec):
    """Preprocess spec to replace all placeholders with actual values"""
    import copy

    # Define all replacements based on spec values
    replacements = {
        'conda_env': spec['conda']['environment_name'],
        'python_version': spec['python']['install_version'],
        'mhkit_python_version': spec['mhkit_python']['version'],
        'version': spec['package']['version']
    }

    # Deep copy spec and replace all placeholders
    processed_spec = copy.deepcopy(spec)
    processed_spec = replace_placeholders_recursively(processed_spec, replacements)

    return processed_spec


def generate_unix_scripts(spec, env, output_dir):
    """Generate Unix (macOS/Linux) shell scripts - both monolithic and step-by-step"""

    # Generate original monolithic script
    template = env.get_template('install_unix.sh.j2')
    content = template.render(**spec)
    output_path = output_dir / "install_mhkit_python_unix.sh"
    with open(output_path, 'w') as f:
        f.write(content)
    output_path.chmod(0o755)
    print(f"Generated: {output_path}")

    # Generate step-by-step scripts
    step_templates = [
        'step1_detect_conda.sh.j2',
        'step2_install_conda.sh.j2',
        'step3_create_env.sh.j2',
        'step4_install_dependencies.sh.j2',
        'step5_post_install.sh.j2'
    ]

    for step_template in step_templates:
        template = env.get_template(step_template)
        content = template.render(**spec)

        # Extract step name from template filename
        step_name = step_template.replace('.sh.j2', '.sh')
        output_path = output_dir / step_name

        with open(output_path, 'w') as f:
            f.write(content)

        # Make executable
        output_path.chmod(0o755)
        print(f"Generated: {output_path}")


def generate_windows_scripts(spec, env, output_dir):
    """Generate Windows PowerShell scripts - both monolithic and step-by-step"""

    # Generate original monolithic script
    template = env.get_template('install_windows.ps1.j2')
    content = template.render(**spec)
    output_path = output_dir / "install_mhkit_python_windows.ps1"
    with open(output_path, 'w') as f:
        f.write(content)
    print(f"Generated: {output_path}")

    # Generate step-by-step scripts
    step_templates = [
        'step1_detect_conda.ps1.j2',
        'step2_install_conda.ps1.j2',
        'step3_create_env.ps1.j2',
        'step4_install_dependencies.ps1.j2',
        'step5_post_install.ps1.j2'
    ]

    for step_template in step_templates:
        template = env.get_template(step_template)
        content = template.render(**spec)

        # Extract step name from template filename
        step_name = step_template.replace('.ps1.j2', '.ps1')
        output_path = output_dir / step_name

        with open(output_path, 'w') as f:
            f.write(content)

        print(f"Generated: {output_path}")


def create_output_directory():
    """Create the shell_scripts output directory"""
    script_dir = Path(__file__).parent
    output_dir = script_dir.parent / "shell_scripts"
    output_dir.mkdir(exist_ok=True)
    return output_dir


def validate_spec(spec):
    """Validate that the spec contains required sections"""
    required_sections = [
        'conda', 'python', 'mhkit_python', 'hooks', 'support'
    ]

    for section in required_sections:
        if section not in spec:
            raise ValueError(f"Missing required section in spec.json: {section}")

    # Validate hooks structure
    hooks = spec['hooks']
    required_hook_sections = ['pre_install', 'post_install', 'environment_setup']
    for hook_section in required_hook_sections:
        if hook_section not in hooks:
            raise ValueError(f"Missing required hook section: {hook_section}")

        if 'commands' not in hooks[hook_section]:
            raise ValueError(f"Missing commands in hook section: {hook_section}")


def main():
    """Main function to generate all scripts"""
    try:
        print("Loading MHKiT specification...")
        spec = load_spec()

        print("Validating specification...")
        validate_spec(spec)

        print("Setting up Jinja2 environment...")
        env = setup_jinja_environment()

        print("Creating output directory...")
        output_dir = create_output_directory()

        print("Preprocessing placeholders...")
        processed_spec = preprocess_spec_for_templates(spec)

        print("Generating shell scripts...")
        generate_unix_scripts(processed_spec, env, output_dir)
        generate_windows_scripts(processed_spec, env, output_dir)

        print("\nScript generation completed successfully!")
        print(f"Unix script: {output_dir}/install_mhkit_python_unix.sh")
        print(f"Windows script: {output_dir}/install_mhkit_python_windows.ps1")
        print(f"Step-by-step scripts also generated in: {output_dir}")

        print("\nThese scripts can now be called from MATLAB's mhkit.install() function.")

    except Exception as e:
        print(f"Error generating scripts: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()