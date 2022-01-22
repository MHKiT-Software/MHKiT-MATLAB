function [data] = cdip_request_parse_workflow(options)
%CDIP_REQUEST_PARSE_WORKFLOW Parses CDIP data from a file or web request
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parses a passed CDIP netCDF file or requests data for a station number
%   (from http://cdip.ucsd.edu/) and parses. Years may be non-consecutive
%   e.g. [2001, 2010]. Time may be sliced by dates (start_date or end date
%   in YYYY-MM-DD). By default 2D variables are not parsed if all 2D
%   variables are needed. See the MHKiT CDIP example live script for
%   information on available parameters.
%   
%   Parameters
%   ----------
%       nc: netCDF variable data
%           netCDF data for the given station number and data type.
%       station_number: string
%           Station number of CDIP wave buoy
%       parameters: string or array of strings
%           Parameters to return. If nan, will return all variables except
%           2D-variables.
%       years: int or array of int
%           Year date, e.g. 2001 or [2001, 2010]
%       start_date: string
%           Start date in YYYY-MM-DD, e.g. '2012-04-01'
%       end_date: string
%           End date in YYYY-MM-DD, e.g. '2012-04-30'
%       data_type: string
%           Either 'historic' or 'realtime'
%       all_2D_variables: boolean
%           Will return all 2D data. Enabling this will add significant
%           processing time. If all 2D variables are not needed it is
%           recommended to pass 2D parameters of interest using the
%           'parameters' keyword and leave this set to False. Default False.
%   
%   Returns
%   -------
%       data: structure
%           'vars1D': structure
%               1D variables indexed by time
%           'metadata': structure
%               Anything not of length time
%           'vars2D': structure of structures, optional
%               If 2D-vars are passed in the 'parameters' key or if run
%               with all_2D_variables=True, then this key will appear
%               with a structure of structures of 2D variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    options.nc = nan;
    options.station_number string = "";
    options.parameters (1,:) string = "";
    options.years (1,:) {mustBeInteger} = nan;
    options.start_date string = "";
    options.end_date string = "";
    options.data_type string {mustBeMember( ...
        options.data_type, {'historic','realtime'})} = "historic";
    options.all_2D_variables logical = false;
end

DATA_CATEGORIES = {'wave', 'sst', 'gps', 'dwr', 'meta'};
DATA_2D_NAMES = {'waveEnergyDensity', 'waveMeanDirection', ...
                 'waveA1Value', 'waveB1Value', ...
                 'waveA2Value', 'waveB2Value', ...
                 'waveCheckFactor', 'waveSpread', ...
                 'waveM2Value', 'waveN2Value'};

if isnan(options.nc) && options.station_number == ""
    throw(MException('MATLAB:cdip_request_parse_workflow', ...
        'Must provide either CDIP netCDF data or a station number.'));
end

% TODO: Should we keep the options.nc parameter? Where would it come from?
if ~isnan(options.nc)
    throw(MException('MATLAB:cdip_request_parse_workflow', ...
        'Optional parameter "nc" is not implemented.'));
end

% Create list of start and end datetimes for which to query data
[start_datetimes, end_datetimes] = start_end_datetimes(options);

% Query info on buoy and available data
%  NOTE: Unlike Python, MATLAB doesn't return all available data from a
%  netCDF query; each variable or attribute needs to be queried separately
url_query = get_url_query(options);
nc_info = ncinfo(url_query);

% Build list of data to query
if options.parameters ~= ""     % if data to query is specified
    data_to_query = options.parameters;
    if options.all_2D_variables == true
        data_to_query = union(data_to_query, DATA_2D_NAMES);   % add all 2D
    end
    % Add 'waveFrequency' if there's any 2D data queried
    if any(ismember(data_to_query, DATA_2D_NAMES))
        data_to_query = union(data_to_query, 'waveFrequency');
    end
    % Add time data for each data category
    categories_in_data = data_categories(data_names, DATA_CATEGORIES);
    categories_in_data = setdiff(categories_in_data, 'meta'); % omit 'meta'
    data_to_query = union(data_to_query, strcat(categories_in_data, 'Time'));
else                            % else query all data except maybe 2D data
    data_to_query = {nc_info.Variables.Name};
    if options.all_2D_variables == false
        data_to_query = setdiff(data_to_query, DATA_2D_NAMES);  % remove 2D
    end
end
data_to_query = sort(data_to_query);

% Compile output structure with queried netCDF data
data.metadata.name = deblank(convertCharsToStrings( ...
    ncread(url_query, 'metaStationName')));
for i = 1:length(start_datetimes)                   % for each time period
    for j = 1:length(data_to_query)                 % for each data metric
        name = data_to_query{j};
        category = data_category(name);
        % TODO: Don't filter each variable as you need to know the
        % associated time. Build structure for each datetime set first.
        % Study Python code again first.
        data_value = get_netcdf_variable( ...
            url_query, name, start_datetimes(i), end_datetimes(i));
        if isnan(data_value)
            warning("MATLAB:cdip_request_parse_workflow", ...
                    "Data name '%s' not found.", name);
        end

        if endsWith(name, '2D')                      % if 2D data
            if ~isfield(data.data.wave2D, name)      % if not yet in struct
                data.data.(strcat(category, '2D')).(name) = data_value;
            else
                % append data
                x=1;
            end
        else                                         % 1D data
            if ~isfield(data.data.wave, name)        % if not yet in struct
                data.data.(category).(name) = data_value;
            else
                % append data
                x=1;
            end
        end
    end
end

data.test = nan;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_netcdf_variable(url_query, data_name, ...
                                    start_datetime, end_datetime)
%GET_NETCDF_VARIABLE Returns variable within start and end time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data_value = ncread(url_query, data_name);
    %TODO: return nan if not found
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function url_query = get_url_query(options)
%GET_URL_QUERY Builds the URL for querying the netCDF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_url = "http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip";
if options.data_type == "historic"
    url_query = sprintf("%s/archive/%sp1/%sp1_historic.nc", ...
                        data_url, ...
                        options.station_number, ...
                        options.station_number);
elseif options.data_type == "realtime"
    url_query = sprintf("%s/realtime/%sp1_rt.nc", ...
                        data_url, ...
                        options.station_number);
end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data_names = data_2D()
% %DATA_2D Returns array of the 2D data names
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     data_names = ['waveEnergyDensity', 'waveMeanDirection', ...
%                   'waveA1Value', 'waveB1Value', ...
%                   'waveA2Value', 'waveB2Value', ...
%                   'waveCheckFactor', 'waveSpread', ...
%                   'waveM2Value', 'waveN2Value'];
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function is_2D = is_data_2D(data_name)
% %IS_DATA_2D Returns whether data is 2D based on its name
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     two_dimensional_data = {'waveEnergyDensity', 'waveMeanDirection', ...
%                             'waveA1Value', 'waveB1Value', ...
%                             'waveA2Value', 'waveB2Value', ...
%                             'waveCheckFactor', 'waveSpread', ...
%                             'waveM2Value', 'waveN2Value'};
%     is_2D = any(ismember(two_dimensional_data, data_name));
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function categories = data_categories(data_names, all_categories)
%DATA_CATEGORIES Returns the category names for the data based on its names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Categories are non-duplicate. Definitely not most efficient algorithm.
wrapper_fun = @(x) data_category(x, all_categories);
category_of_each = cellfun(wrapper_fun, data_names);
categories = unique(category_of_each);
categories = setdiff(categories, '');  % remove any unknown categories ('')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function category = data_category(data_name, all_categories)
%DATA_CATEGORY Returns the category name for the data based on its name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
category = '';
for i = 1:length(all_categories)
    if startsWith(data_name, all_categories{i})     % category is the prefix
        category = all_categories{i};
        break
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [start_datetimes, end_datetimes] = start_end_datetimes(options)
%START_END_DATETIMES Creates list of start and end datetimes to query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_datetimes = {};
end_datetimes = {};
if ~isnan(options.years)
    % Formulate start and end dates from years parameter
    for i = 1:length(options.years)
        start_datetimes{end+1} = datetime(options.years(i), 1, 1, 0, 0, 0, ...
                                          'TimeZone', 'UTC');
        end_datetimes{end+1} = datetime(options.years(i), 12, 31, 23, 59, 59, ...
                                        'TimeZone', 'UTC');
    end
else
    % If start or end date is needed, query times from the netCDF data
    if options.start_date == "" && options.end_date == ""
        waveTime = ncread(url_query, 'waveTime');
    end
    % Substitute in netCDF start/end dates as needed
    if options.start_date ~= ""
        start_datetimes(1) = datetime(options.start_date, ...
                                      'InputFormat', 'yyyy-MM-dd', ...
                                      'TimeZone', 'UTC');
    else
        start_datetimes(1) = datetime(waveTime(1), ...
                                      'ConvertFrom', 'posixtime', ...
                                      'TimeZone', 'UTC');
    end
    if options.end_date ~= ""
        end_datetimes(1) = datetime(options.end_date, ...
                                    'InputFormat', 'yyyy-MM-dd', ...
                                    'TimeZone', 'UTC');
        end_datetimes(1).Hour = 23;
        end_datetimes(1).Minute = 59;
        end_datetimes(1).Second = 59;
    else
        end_datetimes(1) = datetime(waveTime(end), ...
                                    'ConvertFrom', 'posixtime', ...
                                    'TimeZone', 'UTC');
    end
end
end