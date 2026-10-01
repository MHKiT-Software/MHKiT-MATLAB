% Script to create MHKiT-MATLAB toolbox binary, `mhkit_<version>.mltbx`

% Define project properties
project_name = 'mhkit';
% This is unique to this project and first setup in the .prj file. DO NOT CHANGE THIS!
toolbox_guid = '1f426c5c-9e72-4f83-8e42-1e51b296aa29';
display_name = 'Marine and Hydrokinetic Toolkit (MHKiT)';
% MATLAB only supports <Major>.<Minor>.<Patch>.<Build>
% MHKIT-MATLAB uses the following encoding scheme to add version number for non release builds
% Major.Minor.Patch.<Build>
% For release builds version number:
%    * Version should only have Major.Minor.Patch, i.e.: 0.7.0
%
% For dev, alpha, beta and release builds we include an encoded build number
% Where a prefix of
%   9 indicates a dev version
%   8 indicates a alpha version
%   7 indicates a beta version
%   6 indicates a release candidate version
% i.e.: Development build for version 0.7.0-dev.2 would be 0.7.0.9002
project_version = '1.1.0';
matlab_minimum_release_supported = 'R2022b';
summary = 'Marine energy data analysis and visualization toolbox';
description = sprintf([ ...
'MHKiT-MATLAB is a MATLAB package designed for marine energy applications to assist in data processing and visualization. The software package includes functionality for:\n\n' ...
'* Data processing\n' ...
'* Data visualization\n' ...
'* Data quality control\n' ...
'* Resource assessment\n' ...
'* Device performance\n' ...
'* Device loads\n\n' ...
'API documentation: https://mhkit-software.github.io/MHKiT/api.html\n' ...
'Provides comprehensive reference for all functions, classes, and methods in the toolkit.\n\n' ...
'Project documentation: https://mhkit-software.github.io/MHKiT/\n' ...
'Contains detailed guides, examples, and best practices for using MHKiT across all supported languages.\n\n' ...
'Project repository: https://github.com/MHKiT-Software/MHKiT-MATLAB\n' ...
'Project source code, contribute to development, and see the latest updates.\n\n' ...
'Issue tracker: https://github.com/MHKiT-Software/MHKiT-MATLAB/issues\n' ...
'Use this to report bugs, request features, or get help with using the toolkit.\n' ...
]);

author_name = 'Rebecca Fao';
author_email = 'rebecca.fao@nlr.gov';
author_company = 'National Laboratory of the Rockies';

fprintf('Starting %s toolbox build process...\n\n', project_name);

% Set up base paths based on script location
[script_path, ~, ~] = fileparts(mfilename('fullpath'));
project_root = fileparts(script_path);  % Parent directory of this script
fprintf('Project root directory: %s\n', project_root);


fprintf('Setting up paths from %s', project_root);
% Add directories to path
mhkit_path = fullfile(project_root, 'mhkit');
examples_path = fullfile(project_root, 'examples');

% These are the files included in the output MATLAB toolbox
% Unit tests, their data and coverage reports, and the legacy setup.py are excluded
excluded_paths = {fullfile(mhkit_path, 'tests'), fullfile(mhkit_path, 'setup.py')};
is_excluded = @(p) any(cellfun(@(e) strcmp(p, e) || startsWith(p, [e filesep]), excluded_paths));

all_files = dir(fullfile(mhkit_path, '**', '*'));
all_files = all_files(~[all_files.isdir]);
toolbox_files = fullfile({all_files.folder}, {all_files.name});
toolbox_files = toolbox_files(~cellfun(is_excluded, toolbox_files));
toolbox_files = toolbox_files(~endsWith(toolbox_files, {'.DS_Store', '.asv'}));
fprintf('\nIncluding %d files in toolbox\n', numel(toolbox_files));

% Folders added to the user's MATLAB path on install
toolbox_matlab_path = strsplit(genpath(mhkit_path), pathsep);
toolbox_matlab_path = toolbox_matlab_path(~cellfun(@isempty, toolbox_matlab_path));
toolbox_matlab_path = toolbox_matlab_path(~cellfun(is_excluded, toolbox_matlab_path));

fprintf('\nStarting MATLAB toolbox build of %s version %s\n', project_name, project_version);

% Create a new toolbox project
fprintf('Creating %s toolbox project...', display_name);
opts = matlab.addons.toolbox.ToolboxOptions(mhkit_path, toolbox_guid);

% Set basic toolbox properties
fprintf('Configuring toolbox properties...\n');
opts.ToolboxName = display_name;
opts.ToolboxVersion = project_version;
opts.Summary = summary;
opts.Description = description;

% Set author information
opts.AuthorName = author_name;
opts.AuthorEmail = author_email;
opts.AuthorCompany = author_company;

% Set screenshot/logo
fprintf('Checking for toolbox logo...\n');
screenshot_path = fullfile(project_root, 'assets', [project_name '_logo.png']);
if exist(screenshot_path, 'file')
    opts.ToolboxImageFile = screenshot_path;
    fprintf('Logo found and set: %s\n', screenshot_path);
else
    fprintf('Warning: Logo file not found at %s\n', screenshot_path);
end

% Set MATLAB release compatibility
opts.MinimumMatlabRelease = matlab_minimum_release_supported;

% Many modules call MHKiT-Python, which requires a local Python install
opts.SupportedPlatforms.Win64 = true;
opts.SupportedPlatforms.Maci64 = true;
opts.SupportedPlatforms.Glnxa64 = true;
opts.SupportedPlatforms.MatlabOnline = false;

opts.ToolboxFiles = toolbox_files;
opts.ToolboxMatlabPath = toolbox_matlab_path;

% Set output path for the .mltbx file
output_file = fullfile(project_root, sprintf('%s_v%s.mltbx', project_name, project_version));
opts.OutputFile = output_file;
fprintf('\nOutput file will be created at: %s\n', output_file);


% Package the toolbox
fprintf('\nPackaging toolbox...\n');
matlab.addons.toolbox.packageToolbox(opts);
fprintf('\n✓ Toolbox successfully created at: %s\n', output_file);

% Create zip file of examples folder
fprintf('\nCreating zip file of examples folder...\n');
examples_zip_path = fullfile(project_root, sprintf('mhkit_examples_v%s.zip', project_version));
zip(examples_zip_path, examples_path);
fprintf('✓ Examples zip file created at: %s\n', examples_zip_path);

% Create zip file containing pyproject.toml and mhkit_python_utils folder
fprintf('\nCreating zip file of mhkit_python_utils...\n');

% Define paths for Python-related files
python_utils_folder = fullfile(project_root, 'mhkit_python_utils');
pyproject_toml_file = fullfile(project_root, 'pyproject.toml');

% Check if required files/folders exist
if ~exist(python_utils_folder, 'dir')
    error('mhkit_python_utils folder not found at: %s', python_utils_folder);
end

if ~exist(pyproject_toml_file, 'file')
    error('pyproject.toml file not found at: %s', pyproject_toml_file);
end

% Create temporary directory to organize files for zipping
temp_dir = fullfile(project_root, 'temp_python_utils');
if exist(temp_dir, 'dir')
    rmdir(temp_dir, 's');
end
mkdir(temp_dir);

% Copy files to temporary directory
copyfile(python_utils_folder, fullfile(temp_dir, 'mhkit_python_utils'));
copyfile(pyproject_toml_file, fullfile(temp_dir, 'pyproject.toml'));

% Create the zip file
python_utils_zip_path = fullfile(project_root, sprintf('mhkit_python_utils_v%s.zip', project_version));
zip(python_utils_zip_path, temp_dir);

% Clean up temporary directory
rmdir(temp_dir, 's');

fprintf('✓ mhkit_python_utils zip file created at: %s\n', python_utils_zip_path);
