function result = parse_info(env_name, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parse specific fields from 'conda info' command output
%   Extracts active environment, location, Python version, and platform 
%   from conda info output
%
% Parameters
% ------------
%     env_name: string
%         Name of conda environment to get info from
%
%     logger: struct
%         Logger object for error reporting
%
% Returns
% ---------
%     result: struct
%         Structure containing conda environment information:
%         - active_environment: string
%         - active_env_location: string  
%         - python_version: string
%         - platform: string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize result struct
    result = struct();
    
    % 1. Verify the environment exists in conda's environment list
    [status, env_list] = mhkit.sys('conda env list');
    if status ~= 0
        logger.error('Failed to get conda environment list');
        result = struct();
        return;
    end
    
    % Check if our environment name exists in the list
    if ~contains(env_list, env_name)
        logger.error('Environment "%s" not found in conda environment list', env_name);
        result = struct();
        return;
    end
    
    result.active_environment = env_name;
    
    % 2. Get Python version directly from environment
    [status, python_version_raw] = mhkit.sys(sprintf('conda run -n %s python --version', env_name));
    if status ~= 0
        logger.error('Failed to get Python version from environment %s', env_name);
        result = struct();
        return;
    end
    version_match = regexp(python_version_raw, 'Python\s+(\d+\.\d+\.\d+)', 'tokens');
    if isempty(version_match)
        logger.error('Could not parse Python version from: %s', python_version_raw);
        result = struct();
        return;
    end
    result.python_version = version_match{1}{1};
    
    % 3. Get environment location directly from Python
    [status, env_location] = mhkit.sys(sprintf('conda run -n %s python -c "import sys; print(sys.prefix)"', env_name));
    if status ~= 0
        logger.error('Failed to get environment location from environment %s', env_name);
        result = struct();
        return;
    end
    result.active_env_location = strip(env_location);
    
    % 4. Get platform directly from Python
    [status, platform_info] = mhkit.sys(sprintf('conda run -n %s python -c "import platform; print(platform.system().lower() + ''-'' + platform.machine())"', env_name));
    if status ~= 0
        logger.error('Failed to get platform info from environment %s', env_name);
        result = struct();
        return;
    end
    result.platform = strip(platform_info);
    
    % Log success
    logger.info('✓ Retrieved environment info for %s:', env_name);
    logger.info('  Python version: %s', result.python_version);
    logger.info('  Location: %s', result.active_env_location);
    logger.info('  Platform: %s', result.platform);
end
