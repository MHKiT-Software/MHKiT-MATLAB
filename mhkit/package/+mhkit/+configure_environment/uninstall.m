function success = uninstall()
    % Remove MHKiT environment fixes from MATLAB startup script using shared utilities

    % Get spec for constants
    spec = mhkit.spec();

    % Get startup file location using shared utility
    [~, startup_file] = mhkit.matlab.get_matlab_user_directory(spec);

    if ~has_mhkit_fixes(startup_file, spec)
        fprintf('No MHKiT fixes found to remove.\n');
        success = true;
        return;
    end

    success = remove_mhkit_fixes(startup_file, spec);

    if success
        fprintf('✓ Removed MHKiT fixes from: %s\n', startup_file);
        fprintf('Restart MATLAB for changes to take effect.\n');
    else
        fprintf('✗ Failed to remove fixes. Manual edit required: %s\n', startup_file);
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

function success = remove_mhkit_fixes(startup_file, spec)
    % Remove MHKiT fixes using shared utilities

    [read_success, content] = mhkit.io.read_file_safe(startup_file);
    if ~read_success
        success = false;
        return;
    end

    cleaned_content = remove_mhkit_section(content, spec);
    success = mhkit.io.write_file_safe(startup_file, cleaned_content);
end

function cleaned_content = remove_mhkit_section(content, spec)
    % Remove MHKiT section using constants from spec

    start_marker = spec.constants.startup_script.start_marker;
    end_marker = spec.constants.startup_script.end_marker;

    start_pos = strfind(content, start_marker);
    end_pos = strfind(content, end_marker);

    if ~isempty(start_pos) && ~isempty(end_pos) && end_pos > start_pos
        before_section = content(1:start_pos-1);
        after_section = content(end_pos + length(end_marker):end);

        if ~isempty(strtrim(before_section)) && ~isempty(strtrim(after_section))
            cleaned_content = [strtrim(before_section), newline, newline, strtrim(after_section)];
        elseif ~isempty(strtrim(before_section))
            cleaned_content = strtrim(before_section);
        elseif ~isempty(strtrim(after_section))
            cleaned_content = strtrim(after_section);
        else
            cleaned_content = '';
        end
    else
        cleaned_content = content;
    end
end