function dataset = read_NDBC_file(file_name, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Reads a NDBC wave buoy data file (from https://www.ndbc.noaa.gov)
%     into a structure.
%
%     Realtime and historical data files can be loaded with this function.
%
%     Note: With realtime data, missing data is denoted by "MM". With
%     historical data, missing data is denoted using a variable number of
%     9's, depending on the data type (for example: 9999.0 999.0 99.0).
%     'N/A' is automatically converted to missing data.
%
%     Data values are converted to float/int when possible. Column names
%     are also converted to float/int when possible (this is useful when
%     column names are frequency).
%
% Parameters
% ------------
%     file_name : string
%         Name of NDBC wave buoy data file
%
%     missing_value : vector of values (optional)
%         Vector of values that denote missing data. Default is
%         ["MM", 9999, 999, 99] which handles both realtime and historical.
%
% Returns
% ---------
%     dataset : structure
%         For meteorological/standard data:
%             dataset.<ColumnName> : data values for each column
%             dataset.time : datetime values as posixtime (seconds since 1970)
%             dataset.units : structure with units for each column
%
%         For spectral data:
%             dataset.spectrum : spectral density matrix [frequencies x times]
%             dataset.time : datetime values as posixtime (seconds since 1970)
%             dataset.frequency : frequency vector [Hz]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    file_name {mustBeTextScalar, mustBeFile}
end

arguments (Repeating)
    varargin
end

% Parse optional missing values argument
if nargin >= 2
    missing_values = varargin{1};
    if ~isnumeric(missing_values) && ~iscell(missing_values) && ~isstring(missing_values)
        error('MHKiT:read_NDBC_file:InvalidInput', ...
            'missing_value must be a vector of numeric or string values');
    end
else
    % Default missing values for both realtime (MM) and historical (9999, 999, 99)
    missing_values = {9999, 999, 99};
end

% Open file and read header rows
fid = fopen(file_name, 'r');
if fid == -1
    error('MHKiT:read_NDBC_file:FileNotFound', 'Cannot open file: %s', file_name);
end

header_line = fgetl(fid);
units_line = fgetl(fid);

% Parse header - split by whitespace
header = strsplit(strtrim(header_line));

% Check if header is commented (starts with #)
if startsWith(header{1}, '#')
    header{1} = extractAfter(header{1}, '#');
end

% Check if units line exists (second line commented)
units = strsplit(strtrim(units_line));
if startsWith(units{1}, '#')
    units_exist = true;
    units{1} = extractAfter(units{1}, '#');
    % Skip units line - data starts on line 3
    data_start_line = 3;
else
    units_exist = false;
    % No units line - data starts on line 2, rewind to read from line 2
    data_start_line = 2;
end

% Determine date columns - check if minutes column exists
if length(header) >= 5 && strcmpi(header{5}, 'mm')
    date_cols = 5;  % YY MM DD hh mm
else
    date_cols = 4;  % YY MM DD hh
end

% Get data column names (after date columns)
data_header = header(date_cols+1:end);
num_data_cols = length(data_header);

% Get units for data columns (if units exist)
if units_exist
    data_units = units(date_cols+1:end);
end

% Build format string for textscan: date columns + data columns
% Date columns are always integers, data columns are floats (with MM -> NaN)
format_str = repmat('%f ', 1, date_cols + num_data_cols);

% Rewind and skip header lines
frewind(fid);
for i = 1:(data_start_line - 1)
    fgetl(fid);
end

% Read all data using textscan - much faster than readtable
% TreatAsEmpty handles 'MM' and other non-numeric values -> NaN
raw_data = textscan(fid, format_str, 'TreatAsEmpty', {'MM', 'N/A', 'NA'}, ...
    'EmptyValue', NaN, 'CommentStyle', '#');
fclose(fid);

% Extract date columns
year_col = raw_data{1};
month_col = raw_data{2};
day_col = raw_data{3};
hour_col = raw_data{4};

if date_cols == 5
    min_col = raw_data{5};
else
    min_col = zeros(size(year_col));
end

% Convert 2-digit years to 4-digit (assume 1900s for >50, 2000s for <=50)
if max(year_col) < 100
    year_col(year_col > 50) = year_col(year_col > 50) + 1900;
    year_col(year_col <= 50) = year_col(year_col <= 50) + 2000;
end

% Create datetime and convert to posixtime
dt = datetime(year_col, month_col, day_col, hour_col, min_col, 0);
time_posix = posixtime(dt);

% Extract data columns into matrix
num_rows = length(year_col);
data_values = zeros(num_rows, num_data_cols);
for col = 1:num_data_cols
    data_values(:, col) = raw_data{date_cols + col};
end

% Replace missing values with NaN
for i = 1:length(missing_values)
    mv = missing_values{i};
    if isnumeric(mv)
        data_values(data_values == mv) = NaN;
    end
end

% Try to convert column names to float (for spectral/frequency data)
is_spectral = true;
freq_values = zeros(1, num_data_cols);
for i = 1:num_data_cols
    num_val = str2double(data_header{i});
    if isnan(num_val)
        is_spectral = false;
        break;
    end
    freq_values(i) = num_val;
end

% Build output structure
if is_spectral
    % Spectral data format
    % Note: Python returns spectrum as [frequencies x times], so transpose
    dataset.spectrum = data_values';
    dataset.frequency = freq_values';
    dataset.time = time_posix;
else
    % Standard meteorological data format
    % Build struct with specific field order to match original Python wrapper:
    % First data column, units, remaining data columns, time

    % Add first data column
    field_name = matlab.lang.makeValidName(data_header{1});
    dataset.(field_name) = data_values(:, 1);

    % Add units structure if available (as second field)
    if units_exist
        for i = 1:num_data_cols
            field_name_u = matlab.lang.makeValidName(data_header{i});
            if i <= length(data_units)
                dataset.units.(field_name_u) = string(data_units{i});
            end
        end
    end

    % Add remaining data columns
    for i = 2:num_data_cols
        field_name = matlab.lang.makeValidName(data_header{i});
        dataset.(field_name) = data_values(:, i);
    end

    % Add time last
    dataset.time = time_posix;
end

end
