function exists = conda_env_exists(env_name, logger)
    % CONDA_ENV_EXISTS Check if a conda environment exists
    % Returns true if the environment exists, false otherwise

    % Run conda env list and capture the output
    [status, cmdout] = mhkit.sys('conda env list');

    if status ~= 0
        logger.error('Failed to execute conda env list');
    end

    % Check if environment name appears in the output
    % Use regexp to match the exact environment name at start of line
    pattern = ['^' env_name '\s'];
    exists = ~isempty(regexp(cmdout, pattern, 'once', 'lineanchors'));
end
