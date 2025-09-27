function success = configure_startup(spec, logger, auto_configure_mhkit_matlab_python_env)
    % Configure MATLAB startup script with platform-specific environment fixes

    if nargin < 3
        auto_configure_mhkit_matlab_python_env = false;
    end

    success = true;
    platform = mhkit.sys.get_platform();

    % Check if platform has startup fixes
    if ~isfield(spec.hooks.environment_setup, 'startup_fixes') || ...
       ~isfield(spec.hooks.environment_setup.startup_fixes, platform)
        return;
    end

    startup_fixes = spec.hooks.environment_setup.startup_fixes.(platform);
    if isempty(startup_fixes)
        return;
    end

    % Get startup.m location
    [startup_file, matlab_user_dir] = get_startup_file_location();

    % Check if fixes are already applied
    if has_mhkit_fixes(startup_file)
        logger.info('✓ MHKiT startup fixes already configured');
        return;
    end

    % Get user consent (unless auto-configure is enabled)
    if ~auto_configure_mhkit_matlab_python_env
        user_choice = get_user_consent(startup_fixes, startup_file, spec);

        switch user_choice
            case 'n'
                logger.info('Startup script configuration declined - continuing installation');
                return;
            case 'exit'
                success = false;
                return;
        end
    end

    % Apply the fixes
    success = apply_startup_fixes(startup_file, matlab_user_dir, startup_fixes, spec);

    if success
        logger.info('✓ Configured MATLAB startup script: %s', startup_file);
        logger.info('To remove: mhkit.configure_environment.uninstall()');
    else
        logger.warning('Failed to configure startup script - manual setup may be required');
    end
end

function [startup_file, matlab_user_dir] = get_startup_file_location()
    matlab_user_dir = userpath;
    if isempty(matlab_user_dir)
        if ispc
            matlab_user_dir = fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB');
        else
            matlab_user_dir = fullfile(getenv('HOME'), 'Documents', 'MATLAB');
        end
    else
        matlab_user_dir = strtrim(matlab_user_dir);
        if endsWith(matlab_user_dir, pathsep)
            matlab_user_dir = matlab_user_dir(1:end-1);
        end
    end
    startup_file = fullfile(matlab_user_dir, 'startup.m');
end

function has_fixes = has_mhkit_fixes(startup_file)
    has_fixes = false;
    if ~exist(startup_file, 'file'), return; end

    try
        fid = fopen(startup_file, 'r');
        if fid == -1, return; end
        content = fread(fid, '*char')';
        fclose(fid);
        has_fixes = contains(content, '%% MHKiT Environment Setup - START');
    catch
    end
end

function choice = get_user_consent(startup_fixes, startup_file, spec)
    % Process startup fixes for display
    processed_fixes = {};
    for i = 1:length(startup_fixes)
        processed_fixes{i} = substitute_variables(startup_fixes{i}, spec);
    end

    fprintf('\nMHKiT needs to add environment fixes to your MATLAB startup script\n');
    fprintf('to ensure Python libraries work correctly.\n\n');
    fprintf('File: %s\n\n', startup_file);
    fprintf('Code to add:\n');
    for i = 1:length(processed_fixes)
        fprintf('  %s\n', processed_fixes{i});
    end
    fprintf('\nRemove later with: mhkit.configure_environment.uninstall()\n\n');

    while true
        choice = input('Add these fixes? [y/n/exit]: ', 's');
        choice = lower(strtrim(choice));

        if ismember(choice, {'y', 'yes'})
            choice = 'y';
            break;
        elseif ismember(choice, {'n', 'no'})
            choice = 'n';
            break;
        elseif ismember(choice, {'exit', 'quit', 'e', 'q'})
            choice = 'exit';
            break;
        else
            fprintf('Please enter y, n, or exit\n');
        end
    end
end

function success = apply_startup_fixes(startup_file, matlab_user_dir, startup_fixes, spec)
    success = false;
    try
        if ~exist(matlab_user_dir, 'dir')
            mkdir(matlab_user_dir);
        end

        % Read existing content
        existing_content = '';
        if exist(startup_file, 'file')
            fid = fopen(startup_file, 'r');
            if fid ~= -1
                existing_content = fread(fid, '*char')';
                fclose(fid);
            end
        end

        % Remove any existing MHKiT section
        existing_content = remove_existing_mhkit_section(existing_content);

        % Create new MHKiT section
        mhkit_section = create_mhkit_section(startup_fixes, spec);

        % Combine content
        if isempty(existing_content)
            new_content = mhkit_section;
        else
            new_content = [existing_content, newline, newline, mhkit_section];
        end

        % Write startup.m
        fid = fopen(startup_file, 'w');
        if fid == -1, return; end
        fprintf(fid, '%s', new_content);
        fclose(fid);
        success = true;
    catch
    end
end

function cleaned_content = remove_existing_mhkit_section(content)
    start_marker = '%% MHKiT Environment Setup - START';
    end_marker = '%% MHKiT Environment Setup - END';

    start_pos = strfind(content, start_marker);
    end_pos = strfind(content, end_marker);

    if ~isempty(start_pos) && ~isempty(end_pos) && end_pos > start_pos
        before_section = content(1:start_pos-1);
        after_section = content(end_pos + length(end_marker):end);
        cleaned_content = [before_section, after_section];
    else
        cleaned_content = content;
    end
end

function section = create_mhkit_section(startup_fixes, spec)
    section = sprintf('%% MHKiT Environment Setup - START\n');
    section = [section, sprintf('try\n')];

    for i = 1:length(startup_fixes)
        line = startup_fixes{i};
        processed_line = substitute_variables(line, spec);
        section = [section, sprintf('    %s\n', processed_line)];
    end

    section = [section, sprintf('catch\n')];
    section = [section, sprintf('    %% Silently ignore errors\n')];
    section = [section, sprintf('end\n')];
    section = [section, sprintf('%% MHKiT Environment Setup - END\n')];
end

function result = substitute_variables(text, spec)
    result = text;

    if contains(result, '<conda_env>')
        result = strrep(result, '<conda_env>', spec.conda.environment_name);
    end

    if contains(result, '<python_version>')
        result = strrep(result, '<python_version>', spec.python.install_version);
    end

    if contains(result, '<mhkit_python_version>')
        result = strrep(result, '<mhkit_python_version>', spec.mhkit_python.version);
    end

    if contains(result, '<conda_lib_path>')
        conda_env = spec.conda.environment_name;
        [status, conda_info] = mhkit.sys(sprintf('conda run -n %s python -c "import sys; import os; print(os.path.join(os.path.dirname(sys.executable), ''lib''))"', conda_env));
        if status == 0
            conda_lib_path = strtrim(conda_info);
            result = strrep(result, '<conda_lib_path>', conda_lib_path);
        end
    end
end