function success = write_file_safe(filepath, content)
    % Safely write content to file with standardized error handling
    %
    % Parameters:
    %   filepath (string): Path to file to write
    %   content (string): Content to write to file
    %
    % Returns:
    %   success (logical): True if file was written successfully

    success = false;

    try
        % Ensure parent directory exists
        parent_dir = fileparts(filepath);
        if ~isempty(parent_dir) && ~exist(parent_dir, 'dir')
            mkdir(parent_dir);
        end

        fid = fopen(filepath, 'w');
        if fid == -1
            return;
        end

        fprintf(fid, '%s', content);
        fclose(fid);
        success = true;
    catch
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
    end
end