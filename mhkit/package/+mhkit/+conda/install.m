function result = install(commands, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Detects OS and installs miniconda using appropriate commands
%
% Parameters
% ------------
%     commands: struct
%         Structure containing installation commands for different OS:
%         - windows: cell array of commands for Windows
%         - mac.arm: cell array of commands for macOS ARM64
%         - mac.intel: cell array of commands for macOS Intel
%         - linux: cell array of commands for Linux
%
%     logger: struct
%         Logger object for outputting status messages
%
% Returns
% ---------
%     result: logical
%         true if installation succeeds, false if it fails
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get command list based on OS
    if ispc
        % Check admin rights and choose appropriate Windows installation
        has_admin = mhkit.sys.has_admin_rights();
        if has_admin
            logger.info('Administrator rights detected - using admin installation');
            cmd_list = commands.windows.admin;
        else
            logger.info('Limited user rights - using user installation (no PATH modification)');
            cmd_list = commands.windows.user;
        end
    elseif ismac
        if contains(computer('arch'), 'arm64')
            cmd_list = commands.mac.arm;
        else
            cmd_list = commands.mac.intel;
        end
    elseif isunix
        cmd_list = commands.linux;
    else
        logger.error('Unsupported operating system');
        result = false;
    end

    % Execute installation commands with deterministic error handling
    spec = mhkit.spec();

    for cmd = cmd_list
        logger.info('Executing: %s', cmd{1});
        status = system(cmd{1});
        if status ~= 0
            logger.error('Command failed: %s (exit code: %d)', cmd{1}, status);
            logger.error('');
            logger.error('Manual installation required:');
            if ispc
                logger.error('  1. Install conda manually: %s', spec.support.windows_conda_install);
            else
                logger.error('  1. Install conda manually: %s', spec.support.conda_install_instructions);
            end
            logger.error('  2. Report this issue: %s', spec.support.github_issues);
            logger.error('');
            error('Automatic conda installation failed - see links above for manual installation');
        end
    end
    
    % On Windows, add conda to PATH for current session
    if ispc
        % Add common Windows conda paths to MATLAB's PATH
        userprofile = getenv('USERPROFILE');
        conda_paths = {
            fullfile(userprofile, 'miniconda3', 'Scripts')
            fullfile(userprofile, 'miniconda3', 'condabin')
            fullfile(userprofile, 'Miniconda3', 'Scripts') 
            fullfile(userprofile, 'Miniconda3', 'condabin')
            fullfile(userprofile, 'anaconda3', 'Scripts')
            fullfile(userprofile, 'anaconda3', 'condabin')
        };
        
        current_path = getenv('PATH');
        for i = 1:length(conda_paths)
            if exist(conda_paths{i}, 'dir')
                logger.info('Adding to PATH: %s', conda_paths{i});
                new_path = [conda_paths{i} ';' current_path];
                setenv('PATH', new_path);
            end
        end
    end
    
    result = true;
end
