function [matlab_user_dir, startup_file] = get_matlab_user_directory(spec)
    % Get MATLAB user directory and startup.m file path using platform abstraction
    %
    % Parameters:
    %   spec (struct): Specification structure from spec.json
    %
    % Returns:
    %   matlab_user_dir (string): MATLAB user directory path
    %   startup_file (string): Full path to startup.m

    % Try to get from MATLAB first
    matlab_user_dir = userpath;

    if isempty(matlab_user_dir)
        % Fallback to platform-specific defaults from spec
        platform = mhkit.sys.get_platform();

        if strcmp(platform, 'windows')
            template = spec.constants.platform_directories.windows.matlab_user_dir;
            matlab_user_dir = mhkit.io.substitute_env_vars(template);
        else % linux or mac
            template = spec.constants.platform_directories.unix.matlab_user_dir;
            matlab_user_dir = mhkit.io.substitute_env_vars(template);
        end
    else
        % Clean up userpath result
        matlab_user_dir = strtrim(matlab_user_dir);
        if endsWith(matlab_user_dir, pathsep)
            matlab_user_dir = matlab_user_dir(1:end-1);
        end
    end

    startup_file = fullfile(matlab_user_dir, 'startup.m');
end