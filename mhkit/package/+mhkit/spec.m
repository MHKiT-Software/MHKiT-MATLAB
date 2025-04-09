function specs = spec()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reads MHKiT specifications from JSON file
%
% Returns
% ---------
%     specs: struct
%       A comprehensive structure containing software specifications
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the full path to the JSON file
    current_file_path = mfilename('fullpath');
    [current_dir, ~, ~] = fileparts(current_file_path);
    json_path = fullfile(current_dir, 'spec.json');

    % Read the JSON file
    try
        % Read file contents
        file_id = fopen(json_path, 'r');
        raw_json = fread(file_id, inf, 'char=>char')';
        fclose(file_id);

        % Parse JSON
        specs = jsondecode(raw_json);

        % Add method to get directories
        specs.dirs = mhkit.get_dirs();
    catch ME
        % Handle potential errors
        error('Failed to read or parse specifications JSON: %s', ME.message);
    end
end
