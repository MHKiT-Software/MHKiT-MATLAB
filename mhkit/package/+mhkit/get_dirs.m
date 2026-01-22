function dirs = get_mhkit_directories()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get MHKiT application directories for different platforms
%   Returns standardized directory paths for config, data, and cache storage
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     dirs: struct
%         Structure containing directory paths:
%         - config: string (configuration files directory)
%         - data: string (user data directory) 
%         - cache: string (temporary cache directory)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize directories structure
    dirs = struct();

    % Base directory name
    app_name = 'mhkit';

    % Get home directory
    base_home = getenv('HOME');

    % Detect operating system and set directories
    if ispc
        % Windows: Use AppData directories
        dirs.config = fullfile(getenv('APPDATA'), app_name);
        dirs.data = fullfile(getenv('LOCALAPPDATA'), app_name);
        dirs.state = fullfile(getenv('LOCALAPPDATA'), app_name, 'state');
        dirs.cache = fullfile(getenv('LOCALAPPDATA'), app_name, 'cache');

    elseif ismac || isunix
        % XDG-style directory handling for both macOS and Linux

        % Check XDG environment variables
        xdg_config_home = getenv('XDG_CONFIG_HOME');
        xdg_data_home = getenv('XDG_DATA_HOME');
        xdg_state_home = getenv('XDG_STATE_HOME');
        xdg_cache_home = getenv('XDG_CACHE_HOME');

        % Config directory
        if ~isempty(xdg_config_home) && exist(xdg_config_home, 'dir')
            dirs.config = fullfile(xdg_config_home, app_name);
        else
            dirs.config = fullfile(base_home, '.config', app_name);
        end

        % Data directory
        if ~isempty(xdg_data_home) && exist(xdg_data_home, 'dir')
            dirs.data = fullfile(xdg_data_home, app_name);
        else
            dirs.data = fullfile(base_home, '.local', 'share', app_name);
        end

        % State directory
        if ~isempty(xdg_state_home) && exist(xdg_state_home, 'dir')
            dirs.state = fullfile(xdg_state_home, app_name);
        else
            dirs.state = fullfile(base_home, '.local', 'state', app_name);
        end

        % Cache directory
        if ~isempty(xdg_cache_home) && exist(xdg_cache_home, 'dir')
            dirs.cache = fullfile(xdg_cache_home, app_name);
        else
            dirs.cache = fullfile(base_home, '.cache', app_name);
        end

        % Fallback to Application Support if XDG directories don't exist
        if ~exist(dirs.config, 'dir')
            app_support_config = fullfile(base_home, 'Library', 'Application Support', app_name);
            if exist(app_support_config, 'dir')
                dirs.config = app_support_config;
            end
        end

        if ~exist(dirs.data, 'dir')
            app_support_data = fullfile(base_home, 'Library', 'Application Support', app_name);
            if exist(app_support_data, 'dir')
                dirs.data = app_support_data;
            end
        end

        if ~exist(dirs.state, 'dir')
            app_support_state = fullfile(base_home, 'Library', 'Application Support', app_name, 'state');
            if exist(app_support_state, 'dir')
                dirs.state = app_support_state;
            end
        end

        if ~exist(dirs.cache, 'dir')
            app_support_cache = fullfile(base_home, 'Library', 'Caches', app_name);
            if exist(app_support_cache, 'dir')
                dirs.cache = app_support_cache;
            end
        end

    else
        % Fallback for unknown systems
        dirs.config = fullfile(base_home, '.' + app_name + '.config');
        dirs.data = fullfile(base_home, '.' + app_name + '.data');
        dirs.state = fullfile(base_home, '.' + app_name + '.state');
        dirs.cache = fullfile(base_home, '.' + app_name + '.cache');
        warning('Unsupported operating system. Using generic home directory locations.');
    end

    % Ensure all directories exist
    directory_fields = fieldnames(dirs);
    for i = 1:length(directory_fields)
        current_dir = dirs.(directory_fields{i});
        if ~exist(current_dir, 'dir')
            mkdir(current_dir);
        end
    end
end
