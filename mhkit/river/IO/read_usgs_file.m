function datast=read_usgs_file(file_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reads a USGS JSON data file (from https://waterdata.usgs.gov/nwis)
%     into a structure
%
% Parameters
% ----------
%     file_name : str
%         Name of USGS JSON data file
%
% Returns
% -------
%     datast : structure
%
%
%         datast.Data: named according to the parameter's variable description
%
%         datast.time: epoch time [s]
%
%         datast.units: units for each parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    file_name (1,:) char
end
arguments (Output)
    datast struct
end

text = jsondecode(fileread(file_name));

time_series = text.value.timeSeries;
num_series = numel(time_series);

names = strings(num_series, 1);
units = strings(num_series, 1);
series_time = cell(num_series, 1);
series_value = cell(num_series, 1);
keep = false(num_series, 1);

for i = 1:num_series
    try
        series = time_series(i);

        % USGS variableDescription is formatted as "<name>, <units>"
        description = string(series.variable.variableDescription);
        name_parts = split(description, ",");
        names(i) = name_parts(1);
        if numel(name_parts) > 1
            units(i) = name_parts(2);
        end

        % Only the first entry of "values" is used, matching the Python
        % implementation's use of values[0]["value"]
        values = series.values(1).value;

        series_time{i} = local_parse_datetime(string({values.dateTime}'));
        series_value{i} = str2double({values.value}');
        keep(i) = true;
    catch ME
        warning('MATLAB:read_usgs_file', ...
            'Failed to process time series %d: %s', i, ME.message);
    end
end

names = names(keep);
units = units(keep);
series_time = series_time(keep);
series_value = series_value(keep);
num_series = numel(names);

% Union of all timestamps across parameters, mirroring pandas'
% combine_first which outer-joins each series on its datetime index
unique_times = unique(vertcat(series_time{:}));

datast = struct();
for i = 1:num_series
    column = nan(numel(unique_times), 1);
    [is_member, loc] = ismember(series_time{i}, unique_times);
    column(loc(is_member)) = series_value{i}(is_member);

    datast.(names(i)) = column;
    datast.units.(names(i)) = units(i);
end

datast.time = posixtime(unique_times).';

end

function dt = local_parse_datetime(date_strings)
% Parses USGS dateTime strings (ISO 8601, with or without a UTC offset)
% into UTC datetimes, mirroring pandas.to_datetime(..., utc=True).
date_strings = strrep(date_strings, "Z", "+00:00");

has_offset = ~isempty(regexp(date_strings(1), '[+-]\d{2}:\d{2}$', 'once'));

if has_offset
    dt = datetime(date_strings, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSX', ...
        'TimeZone', 'UTC');
else
    dt = datetime(date_strings, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', ...
        'TimeZone', 'UTC');
end
end

