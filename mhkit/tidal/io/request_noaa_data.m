function data=request_noaa_data(station, parameter, start_date, end_date, options)

%     Loads NOAA current data directly from https://api.tidesandcurrents.noaa.gov/api/prod/ using a
%     GET request into a structure. NOAA sets max of 31 days between start and end date.
%     See https://api.tidesandcurrents.noaa.gov/api/prod/ for options. All times are reported as GMT and metric
%     units are returned for data.
%
%     Parameters
%     ----------
%     station : str
%         NOAA current station number (e.g. 's08010')
%
%     parameter : str
%         NOAA paramter (e.g. 'currents')
%
%     start_date : str
%         Start date in the format yyyyMMdd
%
%     end_date : str
%         End date in the format yyyyMMdd
%
%     proxy : None
%         Parameter is now deprecated.
%         To request data from behind a firewall, configure in MATLAB
%         Preferences by navigating to:
%           Home -> Environment -> Preferences
%         then:
%           MATLAB -> Web -> Use a proxy server to connect to the Internet
%         See the following for details:
%           https://www.mathworks.com/help/matlab/import_export/proxy.html
%
%     write_json : str or None (optional)
%         Name of json file to write data
%         to call: request_noaa_data(station,parameter,start_date,end_date,"write_json",write_json)
%
%     Returns
%     -------
%     data : structure
%
%
%         data.id: station ID
%
%         data.name: station name
%
%         data.lat: station Latitude
%
%         data.lon: station Longitude
%
%         data.vars: this will vary depending on parameter input.
%
%         data.time: epoch time [s]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    station string
    parameter string
    start_date string
    end_date string
    options.write_json string = "";
end

MAX_RETRIES = 5;                         % number of query retries if error
REQUIRED_FIELDS = {'id', 'name', 'lat', 'lon', 't'};    % 't' is time
MAX_DAYS_PER_QUERY = 31;

start_datetime = datetime(start_date, ...
        'InputFormat', 'yyyyMMdd', ...
        'TimeZone', 'UTC');
end_datetime = datetime(end_date, ...
        'InputFormat', 'yyyyMMdd', ...
        'TimeZone', 'UTC');

% Query data, splitting period as needed to not exceed NOAA's max days limit
data = struct;      % empty structure
is_first_query = true;
start_period_datetime = start_datetime;
while start_period_datetime <= end_datetime
    end_period_datetime = min(end_datetime, start_period_datetime + days(MAX_DAYS_PER_QUERY - 1));
    end_period = datestr(end_period_datetime, 'yyyymmdd');
    start_period = datestr(start_period_datetime, 'yyyymmdd');

    % Query data in period
    [data_in_period, timeseries_fields] = request_noaa_data_restricted_duration( ...
        station, parameter, start_period, end_period, ...
        MAX_RETRIES, REQUIRED_FIELDS, MAX_DAYS_PER_QUERY);

    if isempty(fieldnames(data_in_period)) || isempty(timeseries_fields)
        % do nothing
    elseif is_first_query
        data = data_in_period;
        is_first_query = false;
    else
        % Append to existing data structure
        for i = 1:length(timeseries_fields)
            data.(timeseries_fields{i}) = cat( ...
                1, ...  % by rows
                data.(timeseries_fields{i}), ...
                data_in_period.(timeseries_fields{i}));
        end
    end

    start_period_datetime = end_period_datetime + days(1);
end

% Write NOAA data to a JSON file
if options.write_json ~= ""
    fid = fopen(options.write_json, 'w');
    fprintf(fid, jsonencode(data));
end

end

function [data, timeseries_fields] = request_noaa_data_restricted_duration( ...
        station, parameter, start_date, end_date, ...
        MAX_RETRIES, REQUIRED_FIELDS, MAX_DAYS_PER_QUERY)
    % Verify dates are less than MAX_DAYS_PER_QUERY apart
    start_datetime = datetime(start_date, ...
        'InputFormat', 'yyyyMMdd', ...
        'TimeZone', 'UTC');
    end_datetime = datetime(end_date, ...
        'InputFormat', 'yyyyMMdd', ...
        'TimeZone', 'UTC');
    msg = 'Date range is greater than %d days in request_noaa_data_restricted_duration()';
    assert(days(end_datetime - start_datetime) + 1 <= MAX_DAYS_PER_QUERY, ...
        sprintf(msg, MAX_DAYS_PER_QUERY));

    % Formulate query
    start_date = convertCharsToStrings(start_date);
    start_date = convertCharsToStrings(start_date);
    start_date = convertCharsToStrings(start_date);
    end_date = convertCharsToStrings(end_date);
    data_url = "https://tidesandcurrents.noaa.gov/api/datagetter?";
    api_query = "begin_date=" + start_date ...
                + "&end_date=" + end_date ...
                + "&station=" + station ...
                + "&product=" + parameter ...
                + "&units=metric&" ...
                + "time_zone=gmt&" ...
                + "application=web_services&" ...
                + "format=xml";

    % Display query
    disp("Data request URL: " + data_url + api_query)

    % Submit query and get data
    for i = 0:MAX_RETRIES
        try
            response = webread(data_url + api_query);
            break;
        catch ME
            if i == MAX_RETRIES
                disp(['MATLAB:request_noaa_data: ', ME.identifier]);
                rethrow(ME)
            else
                pause(1);   % pause(seconds) and retry query
            end
        end
    end

    % Organize a structure containing the metadata and timeseries data
    is_metadata_found = false;
    is_timeseries_keys_found = false;
    data_lines = strsplit(response, '\n');
    for i = 1:length(data_lines)
        if startsWith(data_lines{i}, '<metadata')
            % Parse metadata and initialize output structure using them
            if ~is_metadata_found
                % Note: only key names of one or more letters work
                meta_names = regexp(data_lines{i}, '[a-zA-Z]+(?=\=)', 'match');
                for j = 1:length(meta_names)
                    regex_pattern = strcat('(?<=', meta_names{j}, '\=")(.*?)(?=")');
                    value = regexp(data_lines{i}, regex_pattern, 'match');
                    if meta_names{j} == "lat" || meta_names{j} == "lon"
                        data.(meta_names{j}) = str2double(value{1});
                    else
                        data.(meta_names{j}) = convertCharsToStrings(value{1});
                    end
                end
                is_metadata_found = true;

                % Verify all required metadata is present
                for j = 1:length(REQUIRED_FIELDS)
                    if ~isfield(data, REQUIRED_FIELDS{j})
                        ME = MException('MATLAB:request_noaa_data', ...
                            strcat('Queried data from NOAA is missing field', ...
                                REQUIRED_FIELDS{j}));
                    end
                end
            end
        elseif startsWith(data_lines{i}, '<cu')
            % Parse timeseries data keys and add to output structure
            if ~is_timeseries_keys_found
                % Note: only key names of one or more letters work
                ts_keys = regexp(data_lines{i}, '[a-zA-Z]+(?=\=)', 'match');
                timeseries_fields = [];
                for j = 1:length(ts_keys)
                    if ts_keys{j} == 't'
                        ts_field_struct.(ts_keys{j}) = 'time';
                    else
                        ts_field_struct.(ts_keys{j}) = ts_keys{j};
                    end
                    timeseries_fields = [timeseries_fields, ...
                        convertCharsToStrings(ts_field_struct.(ts_keys{j}))];    % to return
                    data.(ts_field_struct.(ts_keys{j})) = [];
                end
                is_timeseries_keys_found = true;
            end

            % Parse data from current line and add to respective arrays in structure
            for j = 1:length(ts_keys)
                regex_pattern = strcat('(?<=', ts_keys{j}, '\=")(.*?)(?=")');
                value = regexp(data_lines{i}, regex_pattern, 'match');
                if ts_keys{j} == 't'
                    value{1} = posixtime(datetime( ...
                        value{1}, ...
                        'InputFormat', 'yyyy-MM-dd HH:mm', ...
                        'TimeZone', 'UTC'));
                else
                    value{1} = str2double(value{1});
                end
                data.(ts_field_struct.(ts_keys{j}))(end+1, :) = value{1};
            end
        elseif startsWith(data_lines{i}, '<error')
            noaa_error_msg = regexp(data_lines{i}, '(?<=.error.).*(?=..error.)', 'match');
            warning('For station id %s between %s and %s, NOAA returned the error:\n%s', ...
                station, start_date, end_date, noaa_error_msg{1});
            data = struct;      % empty structure
            timeseries_fields = [];
            break
        end
    end
end

