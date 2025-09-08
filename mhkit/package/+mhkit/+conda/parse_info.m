function result = parse_info(env_name, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parse specific fields from 'conda info' command output
%   Extracts active environment, location, Python version, and platform 
%   from conda info output
%
% Parameters
% ------------
%     env_name: string
%         Name of conda environment to get info from
%
%     logger: struct
%         Logger object for error reporting
%
% Returns
% ---------
%     result: struct
%         Structure containing conda environment information:
%         - active_environment: string
%         - active_env_location: string  
%         - python_version: string
%         - platform: string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [status, cmd_out] = mhkit.sys(sprintf('conda run -n %s conda info', env_name));
    if status ~= 0
        logger.error('Failed to execute conda info');
    end

    % Initialize the output struct
    result = struct();

    % Define the fields we want to extract
    fields_to_find = {
        'active environment',
        'active env location',
        'python version',
        'platform'
    };

    disp('Parsing conda info output...');
    disp(cmd_out);

    % Split the input into lines
    lines = regexp(cmd_out, '\n', 'split');

    % Process each line
    for i = 1:length(lines)
        line = strtrim(lines{i});

        % Skip empty lines
        if isempty(line)
            continue;
        end

        % Check if line contains any of our desired fields
        for j = 1:length(fields_to_find)
            field = fields_to_find{j};
            % Extract key part (before the colon)
            parts = regexp(line, ':', 'split');
            if ~isempty(parts)
                line_key = strtrim(parts{1});
                if strcmpi(line_key, field)  % Case-insensitive exact match
                    % Get the value
                    value = strjoin(parts(2:end), ':');  % Rejoin in case value contains colons
                    value = strtrim(value);

                    % Convert field name to valid MATLAB field name
                    field_name = regexprep(lower(field), '\s+', '_');

                    % Store in struct
                    result.(field_name) = value;
                    % fprintf('Matched field "%s" with value "%s"\n', field, value);
                    break;
                end
            end
        end
    end

    % Verify all fields were found
    fprintf('\nFound fields:\n');
    for j = 1:length(fields_to_find)
        field_name = regexprep(lower(fields_to_find{j}), '\s+', '_');
        if isfield(result, field_name)
            fprintf('  %s: %s\n', field_name, result.(field_name));
        else
            fprintf('  WARNING: Field "%s" not found!\n', fields_to_find{j});
        end
    end
end
