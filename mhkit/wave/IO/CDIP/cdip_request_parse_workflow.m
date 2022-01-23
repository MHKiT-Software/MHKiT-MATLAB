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

% Build URL to query
url_query = get_url_query(options);

% Query info on buoy and available data (can't return all vars like Python)
nc_info = ncinfo(url_query);

% Build list of data to query and determine categories
data_to_query = make_data_list(options, nc_info, ...
                               DATA_2D_NAMES, DATA_CATEGORIES);
categories_in_data = data_categories(data_to_query, DATA_CATEGORIES);

% Create list of start and end datetimes for which to query data
datetimes = start_end_datetimes(options);
indices = data_indices(url_query, datetimes, categories_in_data);


% Compile output structure with queried netCDF data
for i = 1:length(data_to_query)                      % for each data metric
    name = data_to_query{i};
    is_2D = any(ismember(DATA_2D_NAMES, name));
    category = data_category(name, DATA_CATEGORIES);
    if is_2D
        category_name = strcat(category, '2D');
    else
        category_name = category;
    end
    N_time_ranges = length(indices.(category).start);
    for j = 1:N_time_ranges
        % Query data
        index_start = indices.(category).start(j);
        index_end = indices.(category).end(j);
        index_count = index_end - index_start + 1;
        % TODO: what about 'meta'?
        % TODO: what if category is ''? Grab all as there are no timestamps?
        try
            if is_2D
                value = ncread(url_query, name, ...
                               [1, index_start], [Inf, index_count]);
            else
                value = ncread(url_query, name, ...
                               index_start, index_count);
            end
        catch
            warning("MATLAB:cdip_request_parse_workflow", ...
                    "Data name '%s' not found.", name);
        end

        % Try adding to existing output field, else create new
        try
            value_in_output = data.data.(category_name).(name);
            data.data.(category_name).(name) = ...
                cat(1, value_in_output, value');     % add rows
        catch
            data.data.(category_name).(name) = value';
        end
    end
end

data.metadata.name = deblank(convertCharsToStrings( ...
    ncread(url_query, 'metaStationName')));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indices = data_indices(url_query, datetime_ranges, categories_in_data)
%DATA_INDICES Returns data indices to query for each category and range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(categories_in_data)
    posixtimes = ncread(url_query, strcat(categories_in_data{i}, 'Time'));
    datetimes = datetime(posixtimes, ...
                         'ConvertFrom', 'posixtime', ...
                         'TimeZone', 'UTC');
    for j = 1:length(datetime_ranges.start)     % for each range
        index_start = find(datetimes>=datetime_ranges.start{j}, 1, 'first');
        index_end = find(datetimes<=datetime_ranges.end{j}, 1, 'last');

        indices.(categories_in_data{i}).start(j) = index_start;
        indices.(categories_in_data{i}).end(j) = index_end;
    end
end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data = get_netcdf_variable(url_query, data_name, ...
%                                     start_datetime, end_datetime)
% %GET_NETCDF_VARIABLE Returns variable within start and end time
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     data_value = ncread(url_query, data_name);
%     %TODO: return nan if not found
% end


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
category_of_each = cellfun(wrapper_fun, data_names, ...
                           'UniformOutput', false);
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
function datetimes = start_end_datetimes(options)
%START_END_DATETIMES Creates structure of start and end datetimes to query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datetimes.start = {};
datetimes.end = {};
if ~isnan(options.years)
    % Formulate start and end dates from years parameter
    for i = 1:length(options.years)
        datetimes.start{end+1} = datetime(options.years(i), 1, 1, 0, 0, 0, ...
                                          'TimeZone', 'UTC');
        datetimes.end{end+1} = datetime(options.years(i), 12, 31, 23, 59, 59, ...
                                        'TimeZone', 'UTC');
    end
else
    % If start or end date is needed, query times from the netCDF data
    if options.start_date == "" && options.end_date == ""
        waveTime = ncread(url_query, 'waveTime');
    end
    % Substitute in netCDF start/end dates as needed
    if options.start_date ~= ""
        datetimes.start(1) = datetime(options.start_date, ...
                                      'InputFormat', 'yyyy-MM-dd', ...
                                      'TimeZone', 'UTC');
    else
        datetimes.start(1) = datetime(waveTime(1), ...
                                      'ConvertFrom', 'posixtime', ...
                                      'TimeZone', 'UTC');
    end
    if options.end_date ~= ""
        datetimes.end(1) = datetime(options.end_date, ...
                                    'InputFormat', 'yyyy-MM-dd', ...
                                    'TimeZone', 'UTC');
        datetimes.end(1).Hour = 23;
        datetimes.end(1).Minute = 59;
        datetimes.end(1).Second = 59;
    else
        datetimes.end(1) = datetime(waveTime(end), ...
                                    'ConvertFrom', 'posixtime', ...
                                    'TimeZone', 'UTC');
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_to_query = make_data_list(options, nc_info, ...
                                        DATA_2D_NAMES, DATA_CATEGORIES)
%MAKE_DATA_LIST Compiles list of data to query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.parameters ~= ""     % if data to query is specified
    data_to_query = options.parameters;
    if options.all_2D_variables == true
        data_to_query = union(data_to_query, DATA_2D_NAMES);   % add all 2D
    end
    % Add 'waveFrequency' if there's any 2D data queried
    if any(ismember(data_to_query, DATA_2D_NAMES))
        data_to_query = union(data_to_query, 'waveFrequency');
    end
    % Add timestamps for each data category
    categories_in_data = data_categories(data_to_query, DATA_CATEGORIES);
    categories_in_data = setdiff(categories_in_data, 'meta'); % omit 'meta'
    data_to_query = union(data_to_query, strcat(categories_in_data, 'Time'));
else                            % else query all data except maybe 2D data
    data_to_query = {nc_info.Variables.Name};
    if options.all_2D_variables == false
        data_to_query = setdiff(data_to_query, DATA_2D_NAMES);  % remove 2D
    end
end
data_to_query = sort(data_to_query);
end
