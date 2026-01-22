function [success, content] = read_file_safe(filepath)
    % Safely read file content with standardized error handling
    %
    % Parameters:
    %   filepath (string): Path to file to read
    %
    % Returns:
    %   success (logical): True if file was read successfully
    %   content (string): File content, empty string on failure

    success = false;
    content = '';

    if ~exist(filepath, 'file')
        return;
    end

    try
        fid = fopen(filepath, 'r');
        if fid == -1
            return;
        end

        content = fread(fid, '*char')';
        fclose(fid);
        success = true;
    catch
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
    end
end