function install()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Install MHKiT Python dependencies and configure MATLAB integration
%
% Parameters
% ------------
%     No parameters required
%
% Returns
% ---------
%     No return values
%         Installs Conda, creates Python environment, installs MHKiT-Python
%         package and utilities, and configures MATLAB-Python integration
%
% Example
% -------
%     mhkit.install()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize logger
    logger = mhkit.utils.get_logger();

    try
        % Log start of installation
        logger.info('Starting installation...');

        spec = mhkit.spec();

        % Step 1: Check MATLAB compatibility
        logger.info('Checking MATLAB compatibility...');
        if mhkit.matlab.less_than(spec.matlab.minimum_version);
            logger.error('Cannot install! MATLAB version older than MHKiT Minimum Supported Version of %s. Please upgrade your MATLAB version to use MHKiT!', spec.matlab.minimum_version);
            return
        end

        if mhkit.matlab.greater_than(spec.matlab.maximum_tested_version);
            logger.warning('MATLAB version newer than MHKiT Newest Supported Version of %s. If stability issues with MHKiT arise please consider downgrading MATLAB!', spec.matlab.maximum_tested_version);
        end

        logger.info('MATLAB version %s is compatiable with MHKiT, moving to next step...', version("-release"));


        % Step 2: Check and install Conda if necessary
        logger.info('Checking Conda installation');
        if ~mhkit.conda.exists();
            logger.info('Installing conda...');
            success = mkhit.conda.install(spec.conda.install, logger);
            if ~success
                logger.error("Failed to install conda");
                return
            end
        else
            logger.info("Found conda, using existing installation");
        end


        % % Step 3: Verify Conda functionality
        % logger.info('Verifying Conda functionality');
        % verify_conda_works(logger);

        % Step 4: Create or verify Conda environment
        conda_env_name = spec.conda.environment_name;
        logger.info('Checking for %s Conda environment', conda_env_name);
        conda_env_exists = mhkit.conda.env_exists(conda_env_name);

        if ~conda_env_exists
            logger.info('Creating %s Conda environment', conda_env_name);
            command = spec.conda.create
            command = replace(command, '<conda_env>', conda_env_name);
            command = replace(command, '<python_version>', spec.python.install_version);
            mhkit.sys(command)
        end


        logger.info('Checking compatability of %s Conda environment', conda_env_name);

        conda_info = mhkit.conda.parse_info(conda_env_name, logger);

        disp(conda_info);

        conda_env_python = conda_info.python_version

        is_conda_python_within_bounds = mhkit.python.version_within(conda_env_python, spec.python.minimum_version, spec.python.maximum_version)

        sprintf("is_conda_python_within_bounds %d\n", is_conda_python_within_bounds);

        if ~is_conda_python_within_bounds
            logger.info('Recreating %s Conda environment', conda_env_name);
            mhkit.sys(sprintf('conda remove -n  %s --all -y', conda_env_name));
            command = spec.conda.create;
            command = replace(command, '<conda_env>', conda_env_name);
            command = replace(command, '<python_version>', spec.python.install_version);
            mhkit.sys(command);
        end

        % Check if `mhkit` is in `conda list`
        has_correct_mhkit_python = false;
        conda_packages = mhkit.conda.list(conda_env_name);

        if isfield(conda_packages, "mhkit")
            if contains(conda_packages.mhkit.version, spec.mhkit_python.version)
                has_correct_mhkit_python = true
                logger.info(sprintf("MHKiT-Python %s already installed...", conda_packages.mhkit.version));
            end
        end


        if ~has_correct_mhkit_python
            logger.info('Using conda to install MHKiT-Python of version %s...', conda_env_name);
            command = spec.mhkit_python.install;
            command = replace(command, '<mhkit_python_version>', spec.mhkit_python.version);
            mhkit.sys(sprintf("conda run -n %s %s", conda_env_name, command));

            % Temporary command to get macos arm to the correct mhkit-python version
            if ismac
                mhkit.sys(sprintf("conda run -n %s pip install --upgrade mhkit==%s", conda_env_name, spec.mhkit_python.version))
            end


            % Step 6: Install or verify MHKiT
            logger.info('Verifying MHKiT-Python version...');
            [status, out] = mhkit.sys(sprintf("conda run -n %s %s", conda_env_name, spec.mhkit_python.verify_version.command));
            mhkit_python_version = strip(out);
            expected_mhkit_python_version = spec.mhkit_python.verify_version.expect;
            if ~contains(mhkit_python_version, expected_mhkit_python_version)
                logger.error("MHKiT-Python version of %s, does not match expected version of %s", mhkit_python_version, expected_mhkit_python_version)
                return
            end

        end

        % Verify MHKiT-Python Operation
        command = spec.mhkit_python.verify_operation.command;
        disp(command);
        command = sprintf("conda run -n %s %s", conda_env_name, command)
        disp(command);
        [status, out] = mhkit.sys(command);
        disp(status);
        disp(out);
        mhkit_python_output = strip(out);
        expected_mhkit_python_output = spec.mhkit_python.verify_operation.expect;
        if ~contains(mhkit_python_output, expected_mhkit_python_output)
            logger.error("MHKiT-Python output of %s, does not match expected output of %s", mhkit_python_output, expected_mhkit_python_output)
            return
        end


        % Step 7: Install or verify mhkit_python_utils
        logger.info('Installing mhkit_python_utils');
        download_path = spec.package.python_package;
        disp(download_path);
        download_path = replace(download_path, "<version>", spec.package.version);
        disp(download_path);
        [status, extracted_path] = mhkit.web.download_and_unzip(download_path, spec.dirs.cache, "mhkit_python_utils")

        disp(status);
        if ~status == 1
            logger.error(sprintf("Failed to download and extract %s", download_path));
        end

        logger.info("Installing mhkit_python_utils")
        mhkit.sys(sprintf("conda run -n %s pip install -e ""%s""", conda_env_name, extracted_path))

        mhkit.sys(sprintf("conda run -n %s python -c ""import mhkit_python_utils; print(mhkit_python_utils.__version__)""", conda_env_name))

        return

        % Step 8: Initialize Python integration
        logger.info('Initializing Python integration');
        initialize_python_integration(environment_name, logger);

        % Final verification
        logger.info('Performing final installation verification');
        verify_installation(logger);

        % Confirm installation with user
        if confirm_installation()
            logger.info('MHKiT installation completed successfully');
        else
            logger.warn('Installation cancelled by user');
            return;
        end

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

function check_matlab_compatibility()
    % Check MATLAB version compatibility
    config = mhkit.spec();
    current_version = version('-release');
    minimum_version = config.matlab.minimum_version;

    % Compare MATLAB versions
    if str2double(current_version(1:4)) < str2double(minimum_version(1:4))
        error('MATLAB version is too low. Minimum required: %s, Current: %s', ...
            minimum_version, current_version);
    end
end

function install_conda_if_needed(logger)
    % Check if Conda is installed, install if not
    try
        % Attempt to run conda command
        [status, ~] = system('conda --version');
        if status ~= 0
            logger.info('Conda not found. Attempting to install.');
            install_conda();
        end
    catch
        logger.info('Conda not found. Attempting to install.');
        install_conda();
    end
end

function install_conda()
    % Placeholder for Conda installation logic
    % This would typically involve downloading and installing Miniconda or Anaconda
    % Actual implementation depends on OS and specific installation method
    error('Conda installation not implemented in this version');
end

function verify_conda_works(logger)
    % Verify Conda is working correctly
    [status, result] = system('conda info');
    if status ~= 0
        logger.error('Conda is not working correctly');
        error('Conda verification failed: %s', result);
    end
end

function env_name = create_or_verify_conda_environment(opts, logger)
    % Create Conda environment if it doesn't exist
    env_name = opts.Environment;

    % Check if environment exists
    [status, ~] = system(sprintf('conda env list | grep %s', env_name));

    if status ~= 0
        % Environment doesn't exist, create it
        logger.info('Creating Conda environment: %s', env_name);
        [status, result] = system(sprintf('conda create -n %s python=%s -y', ...
            env_name, opts.Python));

        if status ~= 0
            logger.error('Failed to create Conda environment');
            error('Environment creation failed: %s', result);
        end
    else
        logger.info('Conda environment %s already exists', env_name);
    end
end

function verify_python_version(env_name, desired_version, logger)
    % Verify Python version in the environment
    cmd = sprintf('conda run -n %s python --version', env_name);
    [status, result] = system(cmd);

    if status ~= 0
        logger.error('Failed to check Python version');
        error('Python version check failed: %s', result);
    end

    % Extract Python version from result
    version_match = regexp(result, 'Python (\d+\.\d+)', 'tokens');
    current_version = version_match{1}{1};

    if ~strcmp(current_version, desired_version)
        logger.warn('Python version mismatch. Desired: %s, Current: %s', ...
            desired_version, current_version);
        % Optionally, recreate environment with correct Python version
    end
end

function install_mhkit(env_name, opts, logger)
    % Install or verify MHKiT
    cmd = sprintf('conda run -n %s pip install mhkit==%s', env_name, opts.Version);
    [status, result] = system(cmd);

    if status ~= 0
        logger.error('Failed to install MHKiT');
        error('MHKiT installation failed: %s', result);
    end
end

function install_mhkit_python_utils(env_name, opts, logger)
    % Install or verify mhkit_python_utils
    cmd = sprintf('conda run -n %s pip install mhkit_python_utils', env_name);
    [status, result] = system(cmd);

    if status ~= 0
        logger.error('Failed to install mhkit_python_utils');
        error('mhkit_python_utils installation failed: %s', result);
    end
end

function initialize_python_integration(env_name, logger)
    % Initialize Python integration
    try
        % Activate the conda environment
        pyenv('Version', sprintf('conda run -n %s which python', env_name));

        % Test Python import
        py.importlib.import_module('mhkit');
        py.importlib.import_module('mhkit_python_utils');

        logger.info('Python integration initialized successfully');
    catch ME
        logger.error('Failed to initialize Python integration');
        error('Python initialization error: %s', ME.message);
    end
end

function verify_installation(logger)
    % Perform final verification of installation
    try
        % Additional checks can be added here
        % For example, running a simple MHKiT function
        py.mhkit.river.performance.tidal_turbine_performance;
        logger.info('Installation verification successful');
    catch ME
        logger.error('Installation verification failed');
        error('Verification error: %s', ME.message);
    end
end

function confirmed = confirm_installation()
    % Prompt user for final confirmation
    response = input('Proceed with installation? (y/n): ', 's');
    confirmed = strcmpi(response, 'y');
end
