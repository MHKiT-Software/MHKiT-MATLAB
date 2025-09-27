function success = uninstall()
    % Remove MHKiT environment fixes from MATLAB startup script

    [startup_file, ~] = get_startup_file_location();

    if ~exist(startup_file, 'file') || ~has_mhkit_fixes(startup_file)
        fprintf('No MHKiT fixes found to remove.\n');
        success = true;
        return;
    end

    success = remove_mhkit_fixes(startup_file);

    if success
        fprintf('✓ Removed MHKiT fixes from: %s\n', startup_file);
        fprintf('Restart MATLAB for changes to take effect.\n');
    else
        fprintf('✗ Failed to remove fixes. Manual edit required: %s\n', startup_file);
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
    try
        fid = fopen(startup_file, 'r');
        if fid == -1, return; end
        content = fread(fid, '*char')';
        fclose(fid);
        has_fixes = contains(content, '%% MHKiT Environment Setup - START');
    catch
    end
end

function success = remove_mhkit_fixes(startup_file)
    success = false;
    try
        fid = fopen(startup_file, 'r');
        if fid == -1, return; end
        content = fread(fid, '*char')';
        fclose(fid);

        cleaned_content = remove_mhkit_section(content);

        fid = fopen(startup_file, 'w');
        if fid == -1, return; end
        fprintf(fid, '%s', cleaned_content);
        fclose(fid);
        success = true;
    catch
    end
end

function cleaned_content = remove_mhkit_section(content)
    start_marker = '%% MHKiT Environment Setup - START';
    end_marker = '%% MHKiT Environment Setup - END';

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