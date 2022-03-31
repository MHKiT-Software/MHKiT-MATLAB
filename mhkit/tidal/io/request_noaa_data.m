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
