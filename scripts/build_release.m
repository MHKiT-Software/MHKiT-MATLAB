% Script to create MHKiT-MATLAB toolbox binary, `mhkit_<version>.mltbx`

% Define project properties
project_name = 'mhkit';
% This is unique to this project and first setup in the .prj file. DO NOT CHANGE THIS!
toolbox_guid = '1f426c5c-9e72-4f83-8e42-1e51b296aa29';
display_name = 'Marine and Hydrokinetic Toolkit (MHKiT)';
project_version = '0.6.0';
matlab_minimum_release_supported = 'R2022b';
summary = 'Marine energy data analysis and visualization toolbox';
description = ['MHKiT-MATLAB is a MATLAB package designed for marine energy applications to assist in\n' ...
    'data processing and visualization.  The software package include functionality for:\n\n' ...
    '* Data processing\n' ...
    '* Data visualization\n' ...
    '* Data quality control\n' ...
    '* Resource assessment\n' ...
    '* Device performance\n' ...
    '* Device loads'];

author_name = 'Rebecca Fao';
author_email = 'rebecca.fao@nrel.gov';
author_company = 'National Renewable Energy Laboratory';

fprintf('Starting %s toolbox build process...\n\n', project_name);

% Set up base paths based on script location
[script_path, ~, ~] = fileparts(mfilename('fullpath'));
project_root = fileparts(script_path);  % Parent directory of this script
fprintf('Project root directory: %s\n', project_root);


fprintf('Setting up paths from %s', project_root);
% Add directories to path
mhkit_path = fullfile(project_root, 'mhkit');
examples_path = fullfile(project_root, 'examples');

% These are the folders and files included in the output MATLAB toolbox
toolbox_files = {mhkitPath};


fprintf('Adding %s to the MATLAB path recursively...\n', mhkitPath)

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
screenshotPath = fullfile(project_root, 'assets', [project_name '_logo.png']);
if exist(screenshotPath, 'file')
    opts.ToolboxImageFile = screenshotPath;
    fprintf('Logo found and set: %s\n', screenshotPath);
else
    fprintf('Warning: Logo file not found at %s\n', screenshotPath);
end

% Set MATLAB release compatibility
opts.MinimumMatlabRelease = matlab_minimum_release_supported;

opts.ToolboxFiles = toolboxFiles;

% Set output path for the .mltbx file
outputFile = fullfile(project_root, sprintf('%s_v%s.mltbx', project_name, project_version));
opts.OutputFile = outputFile;
fprintf('\nOutput file will be created at: %s\n', outputFile);


% Package the toolbox
fprintf('\nPackaging toolbox...\n');
matlab.addons.toolbox.packageToolbox(opts);
fprintf('\n✓ Toolbox successfully created at: %s\n', outputFile);

% Create zip file of examples folder
fprintf('\nCreating zip file of examples folder...\n');
examples_zip_path = fullfile(project_root, sprintf('mhkit_examples_v%s.zip', project_version));
zip(examples_zip_path, examples_path);
fprintf('✓ Examples zip file created at: %s\n', examples_zip_path);
