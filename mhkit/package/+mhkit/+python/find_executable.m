function python_path = find_executable(env_name, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find Python executable in conda environment using multiple methods
%
% Parameters
% ------------
%     env_name: string
%         Name of the conda environment
%     logger: struct
%         Logger object for outputting status messages
%
% Returns
% ---------
%     python_path: string
%         Full path to the Python executable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    logger.info('Locating Python executable for environment "%s"...', env_name);

    % Method 1: conda run with python command
    logger.info('Method 1: Using conda run...');
    python_cmd = mhkit.sys.python_cmd();
    conda_command = sprintf('conda run -n %s %s -c "import sys; print(sys.executable)"', env_name, python_cmd);
    logger.info('Executing: %s', conda_command);
    [status, output] = mhkit.sys(conda_command);
    logger.info('Status: %d, Output: "%s"', status, output);

    if status == 0 && ~isempty(strip(output))
        python_path = strip(output);
        logger.info('✓ Method 1 succeeded: %s', python_path);
        return;
    end

    % Method 2: Direct conda environment path construction
    logger.info('Method 1 failed, trying Method 2: Direct path construction...');
    if ispc
        conda_base = getenv('CONDA_PREFIX');
        if isempty(conda_base)
            conda_base = fullfile(getenv('USERPROFILE'), 'miniconda3');
        end
        python_path = fullfile(conda_base, 'envs', env_name, 'python.exe');
    else
        conda_base = getenv('CONDA_PREFIX');
        if isempty(conda_base)
            conda_base = fullfile(getenv('HOME'), 'miniconda3');
        end
        python_path = fullfile(conda_base, 'envs', env_name, 'bin', 'python');
    end
    logger.info('Constructed path: %s', python_path);

    if exist(python_path, 'file')
        logger.info('✓ Method 2 succeeded: %s', python_path);
        return;
    end

    % Method 3: Parse conda env list output
    logger.info('Method 2 failed, trying Method 3: conda env list...');
    [status, output] = system(sprintf('conda env list | grep %s', env_name));
    logger.info('conda env list status: %d, output: "%s"', status, output);

    if status == 0
        lines = splitlines(strip(output));
        for i = 1:length(lines)
            if contains(lines{i}, env_name)
                parts = strsplit(strip(lines{i}));
                if length(parts) >= 2
                    env_path = parts{end};
                    if ispc
                        python_path = fullfile(env_path, 'python.exe');
                    else
                        python_path = fullfile(env_path, 'bin', 'python');
                    end
                    logger.info('Found environment path: %s', env_path);
                    logger.info('Constructed Python path: %s', python_path);
                    if exist(python_path, 'file')
                        logger.info('✓ Method 3 succeeded: %s', python_path);
                        return;
                    end
                end
            end
        end
    end

    % All methods failed
    error_msg = sprintf('Failed to locate Python executable for conda environment "%s". All discovery methods failed.', env_name);
    mhkit.sys.installation_error(logger, error_msg, 'Python Discovery');
    error(error_msg);
end
