function data=request_noaa_data(station, parameter, start_date, end_date, options)

%     Loads NOAA current data directly from https://tidesandcurrents.noaa.gov/api/ using a 
%     GET request into a structure. NOAA sets max of 31 days between start and end date.
%     See https://co-ops.nos.noaa.gov/api/ for options. All times are reported as GMT and metric
%     units are returned for data.
% 
%     Parameters
%     ----------
%     station : str
%         NOAA current station number (e.g. 's08010')
%     parameter : str
%         NOAA paramter (e.g. 'currents')
%     start_date : str
%         Start date in the format yyyyMMdd
%     end_date : str
%         End date in the format yyyyMMdd 
%     proxy : structure or None (optional)
%         Parameter is now deprecated. Will ignore and print a warning.
%         To request data from behind a firewall, configure in MATLAB Preferences.
%         for example "localhost:8080"
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
    options.proxy string = "";
    options.write_json string = "";
end

if (options.proxy ~= "")
    warning('To use a proxy server, configure in MATLAB Preferences.');
end

REQUIRED_FIELDS = {'id', 'name', 'lat', 'lon', 't'};    % 't' is time

% TODO:
% - split query into 30 day intervals

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
response = webread(data_url + api_query);

% Organize a structure containing the metadata and timeseries data
is_metadata_found = false;
is_keys_found = false;
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
                        strcat('Queried data from NOAA is missing field', REQUIRED_FIELDS{j}));
                end
            end
        end
    elseif startsWith(data_lines{i}, '<cu')
        % Parse data keys and add to output structure
        if ~is_keys_found
            % Note: only key names of one or more letters work
            keys = regexp(data_lines{i}, '[a-zA-Z]+(?=\=)', 'match');
            for j = 1:length(keys)
                if keys{j} == 't'
                    field_names.(keys{j}) = 'time';
                else
                    field_names.(keys{j}) = keys{j};
                end
                data.(field_names.(keys{j})) = [];
            end
            is_keys_found = true;
        end

        % Parse data from current line and add to respective arrays in structure
        for j = 1:length(keys)
            regex_pattern = strcat('(?<=', keys{j}, '\=")(.*?)(?=")');
            value = regexp(data_lines{i}, regex_pattern, 'match');
            if keys{j} == 't'
                value{1} = posixtime(datetime(value{1}, 'InputFormat', 'yyyy-MM-dd HH:mm', 'TimeZone','UTC'));
            else
                value{1} = str2double(value{1});
            end
            data.(field_names.(keys{j}))(end+1, :) = value{1};
        end
    end
end

% Write NOAA data to a JSON file
if options.write_json ~= ""
    fid = fopen(options.write_json, 'w');
    fprintf(fid, jsonencode(data));
end
