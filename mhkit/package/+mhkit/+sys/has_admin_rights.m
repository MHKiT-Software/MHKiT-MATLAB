function has_admin = has_admin_rights()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check if current user has Windows administrator rights
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     has_admin: logical
%         true if user has admin rights, false otherwise
%         Always returns false on non-Windows systems
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    has_admin = false;

    if ~ispc
        return; % Only relevant on Windows
    end

    try
        % Method 1: Try to write to Windows system directory
        test_file = fullfile(getenv('WINDIR'), 'Temp', 'mhkit_admin_test.txt');
        try
            fid = fopen(test_file, 'w');
            if fid ~= -1
                fprintf(fid, 'test');
                fclose(fid);
                delete(test_file);
                has_admin = true;
                return;
            end
        catch
            % Continue to next method
        end

        % Method 2: Check using Windows built-in command
        [status, result] = system('net session >nul 2>&1');
        if status == 0
            has_admin = true;
            return;
        end

        % Method 3: Try using whoami command
        [status, result] = system('whoami /groups | find "S-1-5-32-544" >nul 2>&1');
        if status == 0
            has_admin = true;
            return;
        end

    catch
        % If all methods fail, assume no admin rights (safer default)
        has_admin = false;
    end
end