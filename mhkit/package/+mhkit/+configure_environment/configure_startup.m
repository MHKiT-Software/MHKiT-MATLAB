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

    % Get startup.m location using shared utility
    [matlab_user_dir, startup_file] = mhkit.matlab.get_matlab_user_directory(spec);

    % Check if fixes are already applied
    if has_mhkit_fixes(startup_file, spec)
        logger.info('✓ MHKiT startup fixes already configured');
        return;
    end

    % Get user consent (unless auto-configure is enabled)
    if ~auto_configure_mhkit_matlab_python_env
        user_choice = get_user_consent_for_startup(startup_fixes, startup_file, spec);

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

function has_fixes = has_mhkit_fixes(startup_file, spec)
    % Check if startup.m contains MHKiT fixes using shared utilities

    [success, content] = mhkit.io.read_file_safe(startup_file);
    if ~success
        has_fixes = false;
        return;
    end

    start_marker = spec.constants.startup_script.start_marker;
    has_fixes = contains(content, start_marker);
end

function choice = get_user_consent_for_startup(startup_fixes, startup_file, spec)
    % Get user consent for startup script modification using shared utility

    % Process startup fixes for display
    processed_fixes = {};
    for i = 1:length(startup_fixes)
        processed_fixes{i} = mhkit.io.substitute_variables(startup_fixes{i}, spec);
    end

    % Format message using template from spec
    message = sprintf('%s\n\nFile: %s\n\nCode to add:', ...
        'MHKiT needs to add environment fixes to your MATLAB startup script to ensure Python libraries work correctly.', ...
        startup_file);

    choice = mhkit.ui.get_user_consent(message, processed_fixes);
end

function success = apply_startup_fixes(startup_file, matlab_user_dir, startup_fixes, spec)
    % Apply startup fixes using shared utilities

    % Read existing content
    [read_success, existing_content] = mhkit.io.read_file_safe(startup_file);
    if ~read_success
        existing_content = '';
    end

    % Remove any existing MHKiT section
    cleaned_content = remove_existing_mhkit_section(existing_content, spec);

    % Create new MHKiT section
    mhkit_section = create_mhkit_section(startup_fixes, spec);

    % Combine content
    if isempty(cleaned_content)
        new_content = mhkit_section;
    else
        new_content = [cleaned_content, newline, newline, mhkit_section];
    end

    % Write using shared utility
    success = mhkit.io.write_file_safe(startup_file, new_content);
end

function cleaned_content = remove_existing_mhkit_section(content, spec)
    % Remove existing MHKiT section using constants from spec

    start_marker = spec.constants.startup_script.start_marker;
    end_marker = spec.constants.startup_script.end_marker;

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
    % Create MHKiT section using constants from spec

    start_marker = spec.constants.startup_script.start_marker;
    end_marker = spec.constants.startup_script.end_marker;

    section = sprintf('%s\n', start_marker);
    section = [section, sprintf('try\n')];

    for i = 1:length(startup_fixes)
        line = startup_fixes{i};
        processed_line = mhkit.io.substitute_variables(line, spec);
        section = [section, sprintf('    %s\n', processed_line)];
    end

    section = [section, sprintf('catch\n')];
    section = [section, sprintf('    %% Silently ignore errors\n')];
    section = [section, sprintf('end\n')];
    section = [section, sprintf('%s\n', end_marker)];
end