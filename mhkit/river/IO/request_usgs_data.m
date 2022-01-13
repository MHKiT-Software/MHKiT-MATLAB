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
%         'Instantaneous'.
%         to call: request_usgs_data(station,parameter,start_date.end_date,"data_type",data_type)
%
%     proxy : structure or None (optional)
%         Parameter is now deprecated. Will ignore and print a warning.
%         To request data from behind a firewall, configure in MATLAB Preferences.
%         for example "localhost:8080"
%
%     write_json : str or None (optional)
%         Name of json file to write data
%         to call: request_usgs_data(station,parameter,start_date,end_date,"write_json",write_json)
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
    options.proxy string = "";
    options.write_json string = "";
end

if (options.proxy ~= "")
    warning('To use a proxy server, configure in MATLAB Preferences.');
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
response = webread(data_url + api_query);

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
