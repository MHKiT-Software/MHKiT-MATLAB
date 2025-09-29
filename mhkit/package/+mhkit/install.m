function install(auto_configure_mhkit_matlab_python_env)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Install MHKiT Python dependencies and configure MATLAB integration
%
% Parameters
% ------------
%     auto_configure_mhkit_matlab_python_env (logical, optional):
%         If true, automatically configure startup script without user prompt.
%         Default: false (prompts user for consent)
%
% Returns
% ---------
%     No return values
%         Installs Conda, creates Python environment, installs MHKiT-Python
%         package and utilities, and configures MATLAB-Python integration
%
% Example
% -------
%     mhkit.install()                        % Prompts user for startup config
%     mhkit.install(true)                    % Auto-configures startup script
%     mhkit.install(false)                   % Prompts user for startup config
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Handle optional parameter
    if nargin < 1
        auto_configure_mhkit_matlab_python_env = false;
    end

    % Initialize logger
    logger = mhkit.logging.get_logger();

    try
        fprintf('\nInstalling python dependencies for MHKiT-MATLAB...\n\n');

        spec = mhkit.spec();

        % Step 1: Check MATLAB compatibility
        logger.info('Checking MATLAB compatibility...');
        if mhkit.matlab.less_than(spec.matlab.minimum_version)
            logger.error('Cannot install! MATLAB version older than MHKiT Minimum Supported Version of %s. Please upgrade your MATLAB version to use MHKiT!', spec.matlab.minimum_version);
            return
        end

        if mhkit.matlab.greater_than(spec.matlab.maximum_tested_version)
            logger.warning('MATLAB version newer than MHKiT Newest Supported Version of %s. If stability issues with MHKiT arise please consider downgrading MATLAB!', spec.matlab.maximum_tested_version);
        end

        logger.info('✓ MATLAB version %s is compatible with MHKiT', version("-release"));

        % Step 2: Run platform-specific shell script for Python environment setup
        logger.info('Running Python environment setup script...');
        python_path = run_environment_setup_script(logger);

        if isempty(python_path)
            logger.error('Python environment setup failed');
            return
        end

        logger.info('✓ Python environment setup completed');
        logger.info('Python executable: %s', python_path);

        % Step 3: Install utilities (still done in MATLAB since it requires the package)
        logger.info('Installing utilities...');
        install_utilities(spec, logger);
        logger.info('✓ Utilities installed');

        % Step 4: Configure MATLAB integration
        logger.info('Configuring MATLAB integration...');

        % Execute environment setup hooks for startup fixes
        logger.info('Executing environment setup hooks...');
        spec.auto_configure_mhkit_matlab_python_env = auto_configure_mhkit_matlab_python_env;
        hook_success = mhkit.hooks.execute('environment_setup', spec, logger);
        if ~hook_success
            logger.warning('Environment setup hooks failed, continuing with Python integration');
        else
            logger.info('✓ Environment setup hooks completed');
        end

        % Configure MATLAB's Python integration
        execution_mode = spec.constants.matlab_integration.execution_mode;
        mhkit.python.configure_integration(python_path, execution_mode, logger);

        % Step 5: Final verification
        logger.info('Performing final verification...');
        verify_installation(logger);

        logger.info('\n✓ Installation completed successfully!');
        logger.info('Please restart MATLAB to use MHKiT.\n');

    catch ME
        % Log detailed error information
        logger.error('MHKiT installation failed: %s', ME.message);
        if ~isempty(ME.cause)
            logger.error('Caused by: %s', ME.cause{1}.message);
        end

        % Print stack trace for debugging
        getReport(ME)

        % Rethrow the error
        rethrow(ME);
    end
end


function python_path = run_environment_setup_script(logger)
    % Run step-by-step shell scripts for better user feedback and environment persistence

    % Get the directory containing the shell scripts
    current_file_path = mfilename('fullpath');
    [mhkit_dir, ~, ~] = fileparts(current_file_path);  % This is the +mhkit directory
    [package_dir, ~, ~] = fileparts(mhkit_dir);        % Go up one level to package directory
    scripts_dir = fullfile(package_dir, 'shell_scripts');

    % Determine platform
    if ispc
        logger.info('Installing MHKiT-Python environment (Windows) - Step-by-step...');
        python_path = run_windows_step_by_step_installation(scripts_dir, logger);
    else
        logger.info('Installing MHKiT-Python environment (Unix) - Step-by-step...');
        python_path = run_unix_step_by_step_installation(scripts_dir, logger);
    end
end


function python_path = run_windows_step_by_step_installation(scripts_dir, logger)
    % Run Windows installation using step-by-step scripts for better user feedback

    conda_path = '';

    try
        %% Step 1: Detect Conda
        logger.info('Step 1/5: Detecting conda installation...');
        step1_script = fullfile(scripts_dir, 'step1_detect_conda.ps1');
        if ~exist(step1_script, 'file')
            error('Step 1 script not found: %s', step1_script);
        end

        output = run_powershell_script_with_error_handling(step1_script, 'Step 1', logger);

        conda_detected = parse_script_output(output, 'CONDA_DETECTED', logger);

        if strcmp(conda_detected, 'true')
            conda_path = parse_script_output(output, 'CONDA_PATH', logger);
            logger.info('✓ Found conda at: %s', conda_path);
        else
            logger.info('No conda installation detected, will install miniconda');

            %% Step 2: Install Conda
            logger.info('Step 2/5: Installing miniconda...');
            step2_script = fullfile(scripts_dir, 'step2_install_conda.ps1');
            if ~exist(step2_script, 'file')
                error('Step 2 script not found: %s', step2_script);
            end

            output = run_powershell_script_with_error_handling(step2_script, 'Step 2', logger);

            conda_path = parse_script_output(output, 'CONDA_PATH', logger);
            logger.info('✓ Installed miniconda at: %s', conda_path);
        end

        %% Step 3: Create Environment
        logger.info('Step 3/5: Creating conda environment...');
        step3_script = fullfile(scripts_dir, 'step3_create_env.ps1');
        if ~exist(step3_script, 'file')
            error('Step 3 script not found: %s', step3_script);
        end

        output = run_powershell_script_with_error_handling(step3_script, 'Step 3', logger, conda_path);

        env_exists = parse_script_output(output, 'ENV_EXISTS', logger);
        if strcmp(env_exists, 'false')
            env_created = parse_script_output(output, 'ENV_CREATED', logger);
            if strcmp(env_created, 'true')
                logger.info('✓ Created conda environment: mhkit-matlab-env');
            else
                error('Failed to create conda environment');
            end
        else
            logger.info('✓ Using existing conda environment: mhkit-matlab-env');
        end

        %% Step 4: Install Dependencies
        logger.info('Step 4/5: Installing dependencies (this may take several minutes)...');
        step4_script = fullfile(scripts_dir, 'step4_install_dependencies.ps1');
        if ~exist(step4_script, 'file')
            error('Step 4 script not found: %s', step4_script);
        end

        output = run_powershell_script_with_error_handling(step4_script, 'Step 4', logger, conda_path);

        mhkit_installed = parse_script_output(output, 'MHKIT_INSTALLED', logger);
        utils_installed = parse_script_output(output, 'UTILS_INSTALLED', logger);

        if strcmp(mhkit_installed, 'true') && strcmp(utils_installed, 'true')
            logger.info('✓ Dependencies installed successfully');
        else
            if ~strcmp(mhkit_installed, 'true')
                error('Failed to install mhkit-python');
            end
            if ~strcmp(utils_installed, 'true')
                error('Failed to install mhkit_python_utils');
            end
        end

        %% Step 5: Post-Install Configuration
        logger.info('Step 5/5: Post-install configuration and testing...');
        step5_script = fullfile(scripts_dir, 'step5_post_install.ps1');
        if ~exist(step5_script, 'file')
            error('Step 5 script not found: %s', step5_script);
        end

        output = run_powershell_script_with_error_handling(step5_script, 'Step 5', logger, conda_path);

        functionality_test = parse_script_output(output, 'FUNCTIONALITY_TEST', logger);
        if strcmp(functionality_test, 'passed')
            logger.info('✓ Functionality test passed');
        else
            error('Functionality test failed');
        end

        % Extract Python path
        python_path = parse_script_output(output, 'PYTHON_PATH', logger);
        if isempty(python_path)
            error('Failed to get Python executable path');
        end

        % Verify the Python executable exists
        if ~exist(python_path, 'file')
            error('Python executable not found at reported path: %s', python_path);
        end

        logger.info('✓ Step-by-step installation completed successfully');

    catch ME
        logger.error('Step-by-step installation failed: %s', ME.message);
        python_path = '';
        rethrow(ME);
    end
end


function python_path = run_unix_step_by_step_installation(scripts_dir, logger)
    % Run Unix installation using step-by-step scripts for better user feedback

    conda_path = '';

    try
        %% Step 1: Detect Conda
        logger.info('Step 1/5: Detecting conda installation...');
        step1_script = fullfile(scripts_dir, 'step1_detect_conda.sh');
        if ~exist(step1_script, 'file')
            error('Step 1 script not found: %s', step1_script);
        end

        output = run_script_with_error_handling(step1_script, 'Step 1', logger);

        conda_detected = parse_script_output(output, 'CONDA_DETECTED', logger);

        if strcmp(conda_detected, 'true')
            conda_path = parse_script_output(output, 'CONDA_PATH', logger);
            logger.info('✓ Found conda at: %s', conda_path);
        else
            logger.info('No conda installation detected, will install miniconda');

            %% Step 2: Install Conda
            logger.info('Step 2/5: Installing miniconda...');
            step2_script = fullfile(scripts_dir, 'step2_install_conda.sh');
            if ~exist(step2_script, 'file')
                error('Step 2 script not found: %s', step2_script);
            end

            output = run_script_with_error_handling(step2_script, 'Step 2', logger);

            conda_path = parse_script_output(output, 'CONDA_PATH', logger);
            logger.info('✓ Installed miniconda at: %s', conda_path);
        end

        %% Step 3: Create Environment
        logger.info('Step 3/5: Creating conda environment...');
        step3_script = fullfile(scripts_dir, 'step3_create_env.sh');
        if ~exist(step3_script, 'file')
            error('Step 3 script not found: %s', step3_script);
        end

        output = run_script_with_error_handling(step3_script, 'Step 3', logger, conda_path);

        env_exists = parse_script_output(output, 'ENV_EXISTS', logger);
        if strcmp(env_exists, 'false')
            env_created = parse_script_output(output, 'ENV_CREATED', logger);
            if strcmp(env_created, 'true')
                logger.info('✓ Created conda environment: mhkit-matlab-env');
            else
                error('Failed to create conda environment');
            end
        else
            logger.info('✓ Using existing conda environment: mhkit-matlab-env');
        end

        %% Step 4: Install Dependencies
        logger.info('Step 4/5: Installing dependencies (this may take several minutes)...');
        step4_script = fullfile(scripts_dir, 'step4_install_dependencies.sh');
        if ~exist(step4_script, 'file')
            error('Step 4 script not found: %s', step4_script);
        end

        output = run_script_with_error_handling(step4_script, 'Step 4', logger, conda_path);

        mhkit_installed = parse_script_output(output, 'MHKIT_INSTALLED', logger);
        utils_installed = parse_script_output(output, 'UTILS_INSTALLED', logger);

        if strcmp(mhkit_installed, 'true') && strcmp(utils_installed, 'true')
            logger.info('✓ Dependencies installed successfully');
        else
            if ~strcmp(mhkit_installed, 'true')
                error('Failed to install mhkit-python');
            end
            if ~strcmp(utils_installed, 'true')
                error('Failed to install mhkit_python_utils');
            end
        end

        %% Step 5: Post-Install Configuration
        logger.info('Step 5/5: Post-install configuration and testing...');
        step5_script = fullfile(scripts_dir, 'step5_post_install.sh');
        if ~exist(step5_script, 'file')
            error('Step 5 script not found: %s', step5_script);
        end

        output = run_script_with_error_handling(step5_script, 'Step 5', logger, conda_path);

        functionality_test = parse_script_output(output, 'FUNCTIONALITY_TEST', logger);
        if strcmp(functionality_test, 'passed')
            logger.info('✓ Functionality test passed');
        else
            error('Functionality test failed');
        end

        % Extract Python path
        python_path = parse_script_output(output, 'PYTHON_PATH', logger);
        if isempty(python_path)
            error('Failed to get Python executable path');
        end

        % Verify the Python executable exists
        if ~exist(python_path, 'file')
            error('Python executable not found at reported path: %s', python_path);
        end

        logger.info('✓ Step-by-step installation completed successfully');

    catch ME
        logger.error('Step-by-step installation failed: %s', ME.message);
        python_path = '';
        rethrow(ME);
    end
end


function output = run_script_with_error_handling(script_path, step_name, logger, varargin)
    % DRY helper function to run scripts with consistent error handling and output display

    % Build command with optional arguments
    if nargin > 3
        args = sprintf(' "%s"', varargin{:});
        command = sprintf('bash "%s"%s', script_path, args);
    else
        command = sprintf('bash "%s"', script_path);
    end

    % Run the script
    [status, output] = system(command);

    % Always display script output for debugging
    if ~isempty(output)
        fprintf('\n--- %s Output ---\n%s\n--- End %s Output ---\n', step_name, output, step_name);
    end

    % Handle errors with consistent reporting
    if status ~= 0
        spec = mhkit.spec();
        logger.error('%s failed with exit code: %d', step_name, status);
        logger.error('Script output above. Please report this issue: %s', spec.support.github_issues);
        error('%s failed with exit code: %d', step_name, status);
    end
end


function output = run_powershell_script_with_error_handling(script_path, step_name, logger, varargin)
    % DRY helper function to run PowerShell scripts with consistent error handling and output display

    % Build command with optional arguments
    if nargin > 3
        args = sprintf(' "%s"', varargin{:});
        command = sprintf('powershell -ExecutionPolicy Bypass -File "%s"%s', script_path, args);
    else
        command = sprintf('powershell -ExecutionPolicy Bypass -File "%s"', script_path);
    end

    % Run the script
    [status, output] = system(command);

    % Always display script output for debugging
    if ~isempty(output)
        fprintf('\n--- %s Output ---\n%s\n--- End %s Output ---\n', step_name, output, step_name);
    end

    % Handle errors with consistent reporting
    if status ~= 0
        spec = mhkit.spec();
        logger.error('%s failed with exit code: %d', step_name, status);
        logger.error('Script output above. Please report this issue: %s', spec.support.github_issues);
        error('%s failed with exit code: %d', step_name, status);
    end
end


function value = parse_script_output(output, key, logger)
    % Parse key=value pairs from script output
    value = '';

    % Split output into lines
    lines = strsplit(output, '\n');

    % Look for key=value lines
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if startsWith(line, [key '='])
            value = line((length(key)+2):end);  % Skip 'KEY=' part
            return;
        end
    end

    if isempty(value)
        logger.warning('Could not find %s in script output', key);
    end
end


function python_path = extract_python_path_from_output(output, logger)
    % Extract Python path from shell script output
    % Looking for line like: PYTHON_PATH=/path/to/python

    python_path = '';

    % Split output into lines
    lines = strsplit(output, {'\n', '\r\n'});

    % Look for PYTHON_PATH line
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if startsWith(line, 'PYTHON_PATH=')
            python_path = extractAfter(line, 'PYTHON_PATH=');
            python_path = strtrim(python_path);
            logger.info('Extracted Python path: %s', python_path);
            return;
        end
    end

    logger.warning('Could not find PYTHON_PATH in script output');
end


function install_utilities(spec, logger)
    % Install MHKiT Python utilities package
    % This still needs to be done in MATLAB since it requires downloading the package

    conda_env_name = spec.conda.environment_name;

    % Download utilities package
    download_path = spec.package.python_package;
    download_path = replace(download_path, "<version>", spec.package.version);
    [status, extracted_path] = mhkit.web.download_and_unzip(download_path, spec.dirs.cache, "mhkit_python_utils");

    if ~status == 1
        mhkit.sys.installation_error(logger, sprintf('Failed to download utilities from %s', download_path), 'Utilities Download');
        return
    end

    % Install utilities using pip in the conda environment
    mhkit.sys(sprintf("conda run -n %s pip install -e ""%s""", conda_env_name, extracted_path));

    % Verify installation
    python_cmd = mhkit.sys.python_cmd();
    mhkit.sys(sprintf("conda run -n %s %s -c ""import mhkit_python_utils; print(mhkit_python_utils.__version__)""", conda_env_name, python_cmd));
end


function verify_installation(logger)
    % Perform final verification of installation
    try
        % Simple verification test - check if mhkit modules can be imported
        logger.info('Testing basic mhkit module import...');
        py.importlib.import_module('mhkit');

        % Test a simple function call
        logger.info('Testing basic mhkit functionality...');
        result = py.mhkit.river.performance.circular(30);
        logger.info('Test result: %s', char(result));

        logger.info('Installation verification successful');
    catch ME
        logger.warning('Installation verification failed: %s', ME.message);
        logger.warning('This is common and usually means MATLAB needs to restart to load Python changes');
        logger.info('To test your installation after restarting MATLAB, try: py.mhkit.river.performance.circular(30)');
        logger.info('If you continue to have issues, check that pyenv() points to the correct Python environment');
        % Don't throw error - this is not critical enough to stop installation
    end
end