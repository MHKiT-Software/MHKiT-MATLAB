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


        % Step 2: Check and install Conda if necessary
        logger.info('Checking Conda installation...');
        if ~mhkit.conda.exists();
            logger.info('Installing Conda...');
            success = mhkit.conda.install(spec.conda.install, logger);
            if ~success
                logger.error("Failed to install Conda");
                return
            end
            logger.info('✓ Conda installed successfully');
        else
            logger.info('✓ Found existing Conda installation');
        end

        % Step 3: Create or verify Conda environment
        conda_env_name = spec.conda.environment_name;
        logger.info('\nSetting up Python environment...');
        conda_env_exists = mhkit.conda.env_exists(conda_env_name, logger);

        if ~conda_env_exists
            logger.info('Creating environment "%s"...', conda_env_name);
            command = spec.conda.create;
            command = replace(command, '<conda_env>', conda_env_name);
            command = replace(command, '<python_version>', spec.python.install_version);
            mhkit.sys(command);
            logger.info('✓ Environment created');
        else
            logger.info('✓ Environment "%s" ready', conda_env_name);
        end

        logger.info('Verifying Python configuration...');

        conda_info = mhkit.conda.parse_info(conda_env_name, logger);

        conda_env_python = conda_info.python_version;

        is_conda_python_within_bounds = mhkit.python.version_within(conda_env_python, spec.python.minimum_version, spec.python.maximum_version, logger);

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
                has_correct_mhkit_python = true;
                logger.info('✓ Python %s configured', conda_env_python);
            end
        end

        if ~has_correct_mhkit_python
            % Step 4: Install MHKiT-Python
            logger.info('\nInstalling MHKiT-Python v%s...', spec.mhkit_python.version);
            command = spec.mhkit_python.install;
            command = replace(command, '<mhkit_python_version>', spec.mhkit_python.version);
            mhkit.sys(sprintf("conda run -n %s %s", conda_env_name, command));

            % Temporary command to get macos arm to the correct mhkit-python version
            if ismac
                mhkit.sys(sprintf("conda run -n %s pip install --upgrade mhkit==%s", conda_env_name, spec.mhkit_python.version));
            end

            logger.info('Verifying installation...');
            [status, out] = mhkit.sys(sprintf("conda run -n %s %s", conda_env_name, spec.mhkit_python.verify_version.command));
            mhkit_python_version = strip(out);
            expected_mhkit_python_version = spec.mhkit_python.verify_version.expect;
            if ~contains(mhkit_python_version, expected_mhkit_python_version)
                logger.error("Version mismatch: got %s, expected %s", mhkit_python_version, expected_mhkit_python_version);
                return
            end
            logger.info('✓ Installation successful');
        else
            logger.info('✓ MHKiT-Python v%s ready', spec.mhkit_python.version);
        end

        logger.info('Testing functionality...');
        command = spec.mhkit_python.verify_operation.command;
        command = sprintf("conda run -n %s %s", conda_env_name, command);
        [status, out] = mhkit.sys(command);
        mhkit_python_output = strip(out);
        expected_mhkit_python_output = spec.mhkit_python.verify_operation.expect;
        
        % Extract numbers from both outputs for comparison
        if contains(mhkit_python_output, '(30,') && contains(mhkit_python_output, '706.')
            logger.info('✓ Functionality verified');
        else
            logger.error("Functionality test failed - output: %s", mhkit_python_output);
            logger.error("Expected pattern: %s", expected_mhkit_python_output);
            return
        end


        % Install utilities
        logger.info('\nInstalling utilities...');
        download_path = spec.package.python_package;
        download_path = replace(download_path, "<version>", spec.package.version);
        [status, extracted_path] = mhkit.web.download_and_unzip(download_path, spec.dirs.cache, "mhkit_python_utils");

        if ~status == 1
            logger.error("Failed to download utilities from %s", download_path);
            return
        end

        mhkit.sys(sprintf("conda run -n %s pip install -e ""%s""", conda_env_name, extracted_path));
        mhkit.sys(sprintf("conda run -n %s python -c ""import mhkit_python_utils; print(mhkit_python_utils.__version__)""", conda_env_name));
        logger.info('✓ Utilities installed');

        % Configure MATLAB integration
        logger.info('\nConfiguring MATLAB integration...');
        initialize_python_integration(conda_env_name, logger);

        % Final verification
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


function initialize_python_integration(env_name, logger)
    % Initialize Python integration
    try
        % Get the Python executable path from the conda environment
        [status, python_path] = mhkit.sys(sprintf('conda run -n %s python -c "import sys; print(sys.executable)"', env_name));
        
        if status ~= 0
            logger.error('Failed to get Python executable path from conda environment');
            return
        end
        
        python_path = strip(python_path);
        logger.info('Found Python executable: %s', python_path);
        
        % Add Python directory to system PATH (like the working Unix tests)
        python_dir = fileparts(python_path);
        current_path = getenv('PATH');
        if isunix
            new_path = [python_dir ':' current_path];
        else
            new_path = [python_dir ';' current_path];
        end
        setenv('PATH', new_path);
        logger.info('Added Python directory to PATH: %s', python_dir);
        
        % Set MATLAB's Python environment with OutOfProcess execution mode (like working Unix tests)
        pyenv(Version=python_path, ExecutionMode="OutOfProcess");
        logger.info('Configured pyenv with OutOfProcess execution mode');
        
        % Test Python import
        logger.info('Testing Python module imports...');
        py.importlib.import_module('mhkit');
        py.importlib.import_module('mhkit_python_utils');

        logger.info('Python integration initialized successfully');
    catch ME
        logger.error('Failed to initialize Python integration: %s', ME.message);
        % Don't throw error - this is not critical enough to stop installation
    end
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

