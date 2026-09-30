function data = swan_read_block(swan_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read SAN ASCII block format output and 
% return a struct with modeled values and associated metadata.
%
%     Supports both .mat and .txt file formats.
%
% Parameters
% ------------
%     swan_file : string
%         SWAN file name to import
%
% Returns
% ---------
%     data : structure
%         Structure containing SWAN output variables. Each variable is a
%         sub-structure with fields:
%             .values : matrix of data values [Y x X]
%             .<metadata fields> : metadata from file header (Run, Frame, Unit, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    swan_file {mustBeTextScalar, mustBeFile}
end

% Get file extension
[~, ~, ext] = fileparts(swan_file);
ext = lower(ext);

if strcmp(ext, '.mat')
    % Handle .mat files using native MATLAB load
    data = read_block_mat(swan_file);
else
    % Handle text files (default)
    data = read_block_txt(swan_file);
end

end


function data = read_block_mat(swan_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read SWAN block output saved as a .mat file
%
% Parameters
% ------------
%     swan_file : string
%         SWAN .mat file name to import
%
% Returns
% ---------
%     data : structure
%         One sub-structure per variable in the file, each with a .values
%         matrix [Y x X]. Spaces in variable names become underscores
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mat_data = load(swan_file);

    % Get field names (variables in the .mat file)
    var_names = fieldnames(mat_data);

    % Build output structure with values as matrices
    data = struct();
    for i = 1:length(var_names)
        var_name = var_names{i};
        % Make valid MATLAB field name (replace spaces with underscores)
        safe_name = regexprep(var_name, ' ', '_');
        data.(safe_name).values = mat_data.(var_name);
    end
end


function data = read_block_txt(swan_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read SWAN ASCII block output, one variable block per "% Run" header
%
% Parameters
% ------------
%     swan_file : string
%         SWAN ASCII block file name to import
%
% Returns
% ---------
%     data : structure
%         One sub-structure per variable block, each with a .values matrix
%         [Y x X] and the header metadata fields (Run, Frame, Unit, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fid = fopen(swan_file, 'r');
    if fid == -1
        error('MHKiT:swan_read_block:FileNotFound', 'Cannot open file: %s', swan_file);
    end
    lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    % First pass: find all Run: lines and their positions
    run_lines = [];
    for i = 1:length(lines)
        if startsWith(lines{i}, '% Run')
            run_lines(end+1) = i; %#ok<AGROW>
        end
    end

    % Add end of file marker
    run_lines(end+1) = length(lines) + 1;

    % Initialize output
    data = struct();

    % Process each variable block
    for v = 1:length(run_lines)-1
        start_line = run_lines(v);
        end_line = run_lines(v+1) - 1;

        % Parse this variable block
        [var_name, var_data, var_meta] = parse_variable_block(lines, start_line, end_line);

        if ~isempty(var_name) && ~isempty(var_data)
            % Create safe field name
            safe_name = regexprep(strtrim(var_name), '\s+', '_');
            safe_name = matlab.lang.makeValidName(safe_name);

            % Store data and metadata
            data.(safe_name).values = var_data;
            meta_fields = fieldnames(var_meta);
            for m = 1:length(meta_fields)
                data.(safe_name).(meta_fields{m}) = var_meta.(meta_fields{m});
            end
        end
    end
end


function [var_name, matrix, meta] = parse_variable_block(lines, start_idx, end_idx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parse one variable block of a SWAN ASCII block file into a matrix
% and its header metadata
%
% Parameters
% ------------
%     lines : cell array
%         Every line of the file as a char vector
%     start_idx : double
%         Index of the "% Run" header line that starts the block
%     end_idx : double
%         Index of the last line of the block
%
% Returns
% ---------
%     var_name : char
%         Variable name from the header, empty if the block has no data
%     matrix : matrix
%         Data values [Y x X] with the header unit multiplier applied
%     meta : structure
%         Header metadata fields (Run, Frame, Unit, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    var_name = '';
    matrix = [];
    meta = struct();

    if start_idx > length(lines)
        return;
    end

    % Parse the Run line for metadata
    run_line = lines{start_idx};
    meta = parse_line_metadata(run_line);

    % Get variable name
    if isfield(meta, 'vars')
        var_name = strtrim(meta.vars);
    else
        var_name = sprintf('var_%d', start_idx);
    end

    % Get unit multiplier
    unit_multiplier = 1;
    if isfield(meta, 'Unit')
        unit_str = meta.Unit;
        tokens = regexp(unit_str, '([0-9.Ee+-]+)', 'tokens');
        if ~isempty(tokens)
            unit_multiplier = str2double(tokens{1}{1});
        end
    end
    meta.unitMultiplier = unit_multiplier;

    % Find column header line (after %Y line)
    col_line_idx = 0;
    for i = start_idx+1:min(end_idx, start_idx+10)
        if i <= length(lines) && contains(lines{i}, 'X --->')
            % Column headers are 2 lines after this
            if i+2 <= length(lines)
                col_line_idx = i + 2;
            end
            break;
        end
    end

    % Parse column headers if found
    columns = [];
    if col_line_idx > 0 && col_line_idx <= length(lines)
        col_str = regexprep(lines{col_line_idx}, '^%\s*', '');
        col_parts = strsplit(strtrim(col_str));
        columns = str2double(col_parts);
        columns = columns(~isnan(columns));
    end

    % Find data start (first non-% line after headers)
    data_start = 0;
    for i = start_idx+1:end_idx
        if i <= length(lines) && ~startsWith(lines{i}, '%') && ~isempty(strtrim(lines{i}))
            data_start = i;
            break;
        end
    end

    if data_start == 0
        return;
    end

    % Parse data rows
    row_data = {};
    y_indices = [];

    for i = data_start:end_idx
        if i > length(lines)
            break;
        end
        line = lines{i};

        if startsWith(line, '%')
            continue;
        end

        trimmed = strtrim(line);
        if isempty(trimmed)
            continue;
        end

        % Split on whitespace and periods
        % SWAN block format: "  100 100.100.100.100..."
        % First split by whitespace to get Y index and data string
        parts = regexp(trimmed, '\s+', 'split');
        if length(parts) < 2
            continue;
        end

        y_idx = str2double(parts{1});
        if isnan(y_idx)
            continue;
        end

        % Join remaining parts and split by period
        data_str = strjoin(parts(2:end), '');

        % Split by period (but handle **** for NaN)
        % Replace **** patterns first
        data_str = regexprep(data_str, '\*+', 'NaN.');

        % Split by period
        val_parts = strsplit(data_str, '.');
        val_parts = val_parts(~cellfun(@isempty, val_parts));

        % Convert to numbers
        vals = zeros(1, length(val_parts));
        for k = 1:length(val_parts)
            if strcmpi(val_parts{k}, 'NaN')
                vals(k) = NaN;
            else
                vals(k) = str2double(val_parts{k});
            end
        end

        y_indices(end+1) = y_idx; %#ok<AGROW>
        row_data{end+1} = vals; %#ok<AGROW>
    end

    if isempty(row_data)
        return;
    end

    % Determine matrix size
    n_rows = length(row_data);
    n_cols = max(cellfun(@length, row_data));

    % Build matrix (Y indices go from high to low in file)
    matrix = NaN(n_rows, n_cols);
    for r = 1:n_rows
        vals = row_data{r};
        matrix(r, 1:length(vals)) = vals;
    end

    % Apply unit multiplier
    matrix = matrix * unit_multiplier;
end


function metaDict = parse_line_metadata(line)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parse the key:value pairs of a SWAN "% Run" header line, e.g.
% "% Run:TEST  Frame:  COMPGRID **  Significant wave height, Unit:  0.1000E-01 m"
%
% Parameters
% ------------
%     line : char
%         One "% Run" header line
%
% Returns
% ---------
%     metaDict : structure
%         One field per key (Run, Frame, vars, Unit, etc.) holding the
%         trimmed value text. "**" is read as the key "vars"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    metaDict = struct();

    % Remove leading % and whitespace
    line = regexprep(line, '^%\s*', '');

    % Replace ** with vars:
    line = strrep(line, '**', 'vars:');

    % Replace commas and extra whitespace
    line = regexprep(line, ',', ' ');
    line = regexprep(line, '\s+', ' ');

    % Split by colon
    parts = strsplit(line, ':');

    % Parse key:value pairs
    for i = 1:length(parts)-1
        % Current part ends with key name
        current_words = strsplit(strtrim(parts{i}));
        key = current_words{end};

        % Next part starts with value (everything except last word which is next key)
        next_words = strsplit(strtrim(parts{i+1}));
        if i < length(parts)-1
            % Value is all words except the last (which is the next key)
            val = strjoin(next_words(1:end-1), ' ');
        else
            % Last value - take all words
            val = strjoin(next_words, ' ');
        end

        metaDict.(key) = strtrim(val);
    end
end
