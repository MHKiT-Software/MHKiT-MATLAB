function datast=request_usgs_data(station, parameter, start_date, end_date, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Loads USGS data directly from https://waterdata.usgs.gov/nwis into a structure using a 
%     GET request 
%     
% 
% Parameters
% ----------
%     station : str
%         USGS station number (e.g. '08313000')
%
%     parameter : str
%         USGS paramter ID (e.g. '00060' for Discharge, cubic feet per second)
%
%     start_date : str
%         Start date in the format 'YYYY-MM-DD' (e.g. '2018-01-01')
%
%     end_date : str
%         End date in the format 'YYYY-MM-DD' (e.g. '2018-12-31')
%
%     data_type : str (optional)
%         Data type, options include 'Daily' (return the mean daily value) and 
%         'Instantaneous'. Default = 'Daily'
%         to call: request_usgs_data(station,parameter,start_date.end_date,"data_type",data_type)
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
%         to call: request_usgs_data(station,parameter,start_date,end_date,"write_json",write_json)
%
%     timeout : double (optional)
%         Timeout in seconds for each HTTP request. Default = 30
%         to call: request_usgs_data(station,parameter,start_date,end_date,"timeout",timeout)
%         
% Returns
% -------
%     datast : structure 
%
%
%         datast.Data: named according to the parameter's variable description
%
%         datast.time: datetime
%
%         datast.units: units for each parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    station string
    parameter string
    start_date string
    end_date string
    options.data_type = 'Daily';
    options.write_json string = "";
    options.timeout (1,1) {mustBeNumeric, mustBePositive} = 30;
end

% Formulate query
if strcmp(options.data_type, 'Daily')
    data_url = "https://waterservices.usgs.gov/nwis/dv";
    api_query = "/?format=json&sites=" + station ...
                + "&startDT=" + start_date ...
                + "&endDT=" + end_date ...
                + "&statCd=00003" ...
                + "&parameterCd=" + parameter ...
                + "&siteStatus=all";
else
    data_url = "https://waterservices.usgs.gov/nwis/iv";
    api_query = "/?format=json&sites=" + station ...
                + "&startDT=" + start_date ...
                + "&endDT=" + end_date ...
                + "&parameterCd=" + parameter ...
                + "&siteStatus=all";
end

% Display query
disp("Data request URL: " + data_url + api_query)

% Submit query and get data
web_options = weboptions('Timeout', options.timeout);
response = webread_with_retry(data_url + api_query, web_options);

% Parse data label and units
v = response.value.timeSeries.variable.variableDescription;
vv = strsplit(v, ',');

% Parse data values
vals = {response.value.timeSeries.values.value.value}';
vals = cell2mat(cellfun(@str2num,vals,'uniform',0));

% Parse timestamps corresponding to data
ti = {response.value.timeSeries.values.value.dateTime};

% Organize a structure containing the data label, units, values and timestamps
datast.(vv{1}) = vals;
datast.units.(vv{1}) = string(vv{2});
if strcmp(options.data_type, 'Daily')
    input_format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';
else
    % Note: the trailing 'X' captures the entire 6 character timezone
    %  (e.g., '-08:00') and applies it to the time, so the timezone
    %  should always be 'UTC' for formatting in datetime() or elsewhere
    input_format = 'yyyy-MM-dd''T''HH:mm:ss.SSSX';
end
datast.time = posixtime( ...
    datetime( ...
        ti, ...
        'TimeZone', 'UTC', ...
        'InputFormat', input_format));

% Write USGS data to a JSON file
if options.write_json ~= ""
    fid = fopen(options.write_json, 'w');
    fprintf(fid, jsonencode(response));
end

end


function response = webread_with_retry(url, web_options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run a USGS web request, retrying transient failures with exponential
% backoff and jitter
%
% Server errors (HTTP 5xx), rate limiting (HTTP 429), request timeouts
% (HTTP 408), and connection failures are transient, so the request is
% retried after 1, 2, 4, 8, 16, and 32 seconds, each scaled by a random
% factor between 0.5 and 1.5 so repeated clients do not retry in
% lockstep. Other client errors (HTTP 4xx, e.g. an invalid station or
% parameter) are not transient and are raised immediately.
%
% Parameters
% ------------
%     url : string
%         USGS request URL
%     web_options : weboptions
%         Options passed to webread
%
% Returns
% ---------
%     response : struct
%         Decoded JSON response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_RETRIES = 6;
BASE_DELAY = 1;   % seconds, doubled after each failed attempt

for attempt = 0:MAX_RETRIES
    try
        response = webread(url, web_options);
        return
    catch ME
        status = regexp(ME.identifier, 'HTTP(\d+)StatusCodeError', 'tokens', 'once');
        if isempty(status)
            is_transient = true;   % timeout or connection failure
        else
            status = str2double(status{1});
            is_transient = status >= 500 || status == 429 || status == 408;
        end

        if ~is_transient || attempt == MAX_RETRIES
            disp(['MATLAB:request_usgs_data: ', ME.identifier]);
            rethrow(ME)
        end

        delay = BASE_DELAY * 2^attempt * (0.5 + rand());
        fprintf('USGS request failed (%s), retry %d of %d in %.1f s\n', ...
            ME.identifier, attempt + 1, MAX_RETRIES, delay);
        pause(delay);
    end
end
end
