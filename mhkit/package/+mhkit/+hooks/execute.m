function success = execute(hook_name, spec, logger)
    % Execute platform-specific hooks from spec.json
    %
    % Parameters:
    %   hook_name (string): Name of the hook to execute ('pre_install', 'post_install', 'environment_setup')
    %   spec (struct): Specification structure from spec.json
    %   logger (struct): Logger instance
    %
    % Returns:
    %   success (logical): True if all hook commands succeeded
    
    success = true;
    
    % Check if hooks section exists
    if ~isfield(spec, 'hooks')
        logger.debug('No hooks section found in spec');
        return;
    end
    
    % Check if specific hook exists
    if ~isfield(spec.hooks, hook_name)
        logger.debug('Hook "%s" not found in spec.hooks', hook_name);
        return;
    end
    
    % Determine platform
    platform = mhkit.sys.get_platform();
    hook_section = spec.hooks.(hook_name);
    
    logger.info('Executing %s hooks for platform "%s"...', hook_name, platform);
    
    % Execute environment variables first
    if isfield(hook_section, 'environment_variables')
        success = execute_environment_variables(hook_section.environment_variables, platform, spec, logger);
        if ~success
            logger.error('Environment variable setup failed');
            return;
        end
    end
    
    % Execute commands
    if isfield(hook_section, 'commands')
        success = execute_commands(hook_section.commands, platform, spec, logger);
        if ~success
            logger.error('Command execution failed');
            return;
        end
    end

    % Execute startup fixes (for environment_setup hook only)
    if strcmp(hook_name, 'environment_setup') && isfield(hook_section, 'startup_fixes')
        success = execute_startup_fixes(hook_section.startup_fixes, platform, spec, logger);
        if ~success
            logger.error('Startup script configuration failed');
            return;
        end
    end

    logger.info('Successfully completed %s hooks for platform "%s"', hook_name, platform);
end

function success = execute_environment_variables(env_vars_section, platform, spec, logger)
    % Execute environment variable setup for a platform
    %
    % Parameters:
    %   env_vars_section (struct): Environment variables section from spec
    %   platform (string): Target platform
    %   spec (struct): Specification structure
    %   logger (struct): Logger instance
    %
    % Returns:
    %   success (logical): True if all environment variables were set
    
    success = true;
    
    % Check if platform-specific environment variables exist
    if ~isfield(env_vars_section, platform)
        logger.debug('No environment variables found for platform "%s"', platform);
        return;
    end
    
    platform_env_vars = env_vars_section.(platform);
    
    % Skip if empty
    if isempty(platform_env_vars) || (isstruct(platform_env_vars) && isempty(fieldnames(platform_env_vars)))
        logger.debug('No environment variables to set for platform "%s"', platform);
        return;
    end
    
    % Set each environment variable
    var_names = fieldnames(platform_env_vars);
    for i = 1:length(var_names)
        var_name = var_names{i};
        var_value = platform_env_vars.(var_name);
        
        % Perform variable substitution
        processed_value = mhkit.io.substitute_variables(var_value, spec);
        
        logger.debug('Setting environment variable: %s=%s', var_name, processed_value);
        
        try
            setenv(var_name, processed_value);
            logger.info('Set environment variable: %s=%s', var_name, processed_value);
        catch ME
            logger.error('Failed to set environment variable %s: %s', var_name, ME.message);
            success = false;
            return;
        end
    end
end

function success = execute_commands(commands_section, platform, spec, logger)
    % Execute shell commands for a platform
    %
    % Parameters:
    %   commands_section (struct): Commands section from spec
    %   platform (string): Target platform
    %   spec (struct): Specification structure
    %   logger (struct): Logger instance
    %
    % Returns:
    %   success (logical): True if all commands succeeded
    
    success = true;
    
    % Check if platform-specific commands exist
    if ~isfield(commands_section, platform)
        logger.debug('No commands found for platform "%s"', platform);
        return;
    end
    
    platform_commands = commands_section.(platform);
    
    % Skip if empty
    if isempty(platform_commands)
        logger.debug('No commands to execute for platform "%s"', platform);
        return;
    end
    
    % Execute each command
    for i = 1:length(platform_commands)
        command = platform_commands{i};
        success = execute_single_command(command, spec, logger);
        if ~success
            logger.error('Command failed: %s', command);
            return;
        end
    end
end

function success = execute_single_command(command, spec, logger)
    % Execute a single command with variable substitution
    %
    % Parameters:
    %   command (string): Command to execute
    %   spec (struct): Specification structure for variable substitution
    %   logger (struct): Logger instance
    %
    % Returns:
    %   success (logical): True if command succeeded
    
    % Perform variable substitution
    processed_command = mhkit.io.substitute_variables(command, spec);
    
    logger.debug('Executing command: %s', processed_command);
    
    % Execute command
    [status, result] = mhkit.sys(processed_command);
    
    if status == 0
        logger.debug('Command succeeded: %s', processed_command);
        success = true;
    else
        logger.error('Command failed with status %d: %s', status, processed_command);
        if ~isempty(result)
            logger.error('Error output: %s', result);
        end
        success = false;
    end
end

function success = execute_startup_fixes(startup_fixes_section, platform, spec, logger)
    % Execute startup script configuration for a platform

    success = true;

    % Check if platform-specific startup fixes exist
    if ~isfield(startup_fixes_section, platform)
        logger.debug('No startup fixes found for platform "%s"', platform);
        return;
    end

    platform_startup_fixes = startup_fixes_section.(platform);

    % Skip if empty
    if isempty(platform_startup_fixes)
        logger.debug('No startup fixes to apply for platform "%s"', platform);
        return;
    end

    % Call the modular configure_environment function
    % Get auto_configure setting from spec (default to false if not present)
    auto_configure = false;
    if isfield(spec, 'auto_configure_mhkit_matlab_python_env')
        auto_configure = spec.auto_configure_mhkit_matlab_python_env;
    end

    success = mhkit.configure_environment.configure_startup(spec, logger, auto_configure);
end