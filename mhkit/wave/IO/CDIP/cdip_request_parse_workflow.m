function [data] = cdip_request_parse_workflow(options)

%CDIP_REQUEST_PARSE_WORKFLOW Parses CDIP data from a web request
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Requests data for a station number (from http://cdip.ucsd.edu/) and
%   parses. Years may be non-consecutive e.g. [2001, 2010]. Time may be
%   sliced by dates (start_date or end date in YYYY-MM-DD). By default 2D
%   variables are not parsed if all 2D variables are needed.
%
%   Parameters
%   ----------
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
%           'parameters' keyword and leave this set to the default false.
%
%   Returns
%   -------
%       data: structure
%           'data': structure
%               grouped 1D and 2D structures of array data with datetimes
%           'metadata': structure
%               Anything not of length time, including buoy name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    options.station_number string;
    options.parameters (1,:) string = "";
    options.years (1,:) {mustBeInteger} = -1;
    options.start_date string = "";
    options.end_date string = "";
    options.data_type string {mustBeMember( ...
        options.data_type, {'historic','realtime'})} = "historic";
    options.all_2D_variables logical = false;
end

DATA_GROUPS = {'wave', 'sst', 'gps', 'dwr', 'meta'};

% Build URL to query
url_query = get_url_query(options);

% Query info on buoy and available data (can't return all vars like Python)
% converted to table to query like: nc_info.Variables{'waveTime', 'Size'}{1};
nc_info = ncinfo_autoretry(url_query);
nc_info.Variables = struct2table(nc_info.Variables);
nc_info.Variables.Properties.RowNames = nc_info.Variables.Name;

% Build list of data to query
data_to_query = make_data_list(options, nc_info, DATA_GROUPS);

% Create list of start and end datetimes/indices for which to query data
datetimes = start_end_datetimes(options);
indices = data_indices(url_query, datetimes, data_to_query, DATA_GROUPS);

% Query data and compile into output structure
for i = 1:length(data_to_query)                     % for each data metric
    name = data_to_query{i};
    [type, group, shape] = datum_categories(name, nc_info, DATA_GROUPS);
    if shape == "2D"
        group_name = strcat(group, '2D');
    else
        group_name = group;
    end

    if type ~= "data" || shape == "0D"
        % Query it all and add to output
        try
            value = ncread_autoretry(url_query, name);
            data.(type).(group_name).(name) = value;
        catch
            warning("MATLAB:cdip_request_parse_workflow", ...
                    "Data name '%s' not found.", name);
        end
    elseif type == "data" && (shape == "2D" || shape == "1D")
        % number of group time ranges within desired time ranges, if any
        if isfield(indices, group)
            N_time_ranges = length(indices.(group).start);
        else
            N_time_ranges = 0;
        end

        for j = 1:N_time_ranges
            % Query data
            index_start = indices.(group).start(j);
            index_end = indices.(group).end(j);
            index_count = index_end - index_start + 1;
            try
                if shape == "2D"
                    value = ncread_autoretry(url_query, name, ...
                                   [1, index_start], [Inf, index_count]);
                    value = value';
                elseif shape == "1D"
                    value = ncread_autoretry(url_query, name, ...
                                   index_start, index_count);
                end
            catch ME
                if ME.identifier == "MATLAB:imagessci:netcdf:libraryFailure"
                    warning("MATLAB:cdip_request_parse_workflow", ...
                            "Access failure to NetCDF file when querying %s.", name);
                else
                    warning("MATLAB:cdip_request_parse_workflow", ...
                            "Data name '%s' not found.", name);
                end
                continue;
            end

            % Convert any times
            if endsWith(name, 'Time')
                value = datetime(value, ...
                                 'ConvertFrom', 'posixtime', ...
                                 'TimeZone', 'UTC');
            end

            % Try adding to existing output field, else create new
            try
                value_in_output = data.(type).(group_name).(name);
                data.data.(group_name).(name) = ...
                    cat(1, value_in_output, value);     % add rows
            catch
                data.(type).(group_name).(name) = value;
            end
        end
    end
end

% Add buoy name to output
data.metadata.name = deblank(convertCharsToStrings( ...
    ncread_autoretry(url_query, 'metaStationName')));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function group = data_group(data_name, all_groups)
%DATA_GROUP Returns the group name for the data based on its name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group = 'other';
for i = 1:length(all_groups)
    if startsWith(data_name, all_groups{i})     % group is the prefix
        group = all_groups{i};
        break
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function groups = data_groups(data_names, all_groups)
%DATA_GROUPS Returns the group set in data_names based on its names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groups are non-duplicate. Definitely not most efficient algorithm.
wrapper_fun = @(x) data_group(x, all_groups);
group_of_each = cellfun(wrapper_fun, data_names, ...
                           'UniformOutput', false);
groups = unique(group_of_each);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function indices = data_indices(url_query, datetime_ranges, ...
                                data_to_query, all_groups)
%DATA_INDICES Returns data indices to query for each group and range
% e.g.,  indices.wave.start = <index>
%        indices.wave.end = <index>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groups_in_data = data_groups(data_to_query, all_groups);

indices = struct;
for i = 1:length(groups_in_data)
    posixtimes = ncread_autoretry( ...
            url_query, strcat(groups_in_data{i}, 'Time'));
    if isscalar(posixtimes) && isnan(posixtimes)
        continue
    end

    datetimes = datetime(posixtimes, ...
                         'ConvertFrom', 'posixtime', ...
                         'TimeZone', 'UTC');
    for j = 1:length(datetime_ranges.start)     % for each range
        index_start = find(datetimes>=datetime_ranges.start{j}, 1, 'first');
        index_end = find(datetimes<=datetime_ranges.end{j}, 1, 'last');

        if ~isempty(index_start) && ~isempty(index_end)
            indices.(groups_in_data{i}).start(j) = index_start;
            indices.(groups_in_data{i}).end(j) = index_end;
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [type, group, shape] = datum_categories(datum_name, nc_info, all_groups)
%DATUM_CATEGORIES Returns the categorizations for the datum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% type = {'data', 'metadata'}
% group = {'wave', 'sst', 'gps', 'dwr', 'other'}
% shape = {'0D', '1D', '2D'}

% Determine group
group = data_group(datum_name, all_groups);

% Determine shape
size = nc_info.Variables{datum_name, 'Size'}{1};
datatype = nc_info.Variables{datum_name, 'Datatype'}{1};
if length(size) == 2
    shape = '2D';
elseif length(size) == 1 && size(1) > 1 && ...
        datatype ~= "char" && datatype ~= "string"
    shape = '1D';
else
    shape = '0D';
end

% Determine type
try
    time_length = nc_info.Variables{strcat(group, 'Time'), 'Size'}{1};
    if group ~= "other" && size(end) == time_length
        type = 'data';
    else
        type = 'metadata';
    end
catch
    type = 'metadata';      % if there is no '..Time' metric in group
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_2D_names = find_data_2D_names(nc_info)
%FIND_DATA_2D_NAMES Finds names of available 2D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_freq = nc_info.Variables{'waveFrequency', 'Size'}{1};
N_time = nc_info.Variables{'waveTime', 'Size'}{1};
data_2D_names = {};
for i = 1:height(nc_info.Variables)
    if isequal(nc_info.Variables.Size{i}, [N_freq, N_time])
        data_2D_names{end+1} = nc_info.Variables.Name{i};
    end
end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_to_query = make_data_list(options, nc_info, DATA_GROUPS)
%MAKE_DATA_LIST Compiles list of data to query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_2D_names = find_data_2D_names(nc_info);

if options.parameters ~= ""               % if data to query is specified
    % Exclude non-existent parameters
    data_to_query = intersect( ...
        options.parameters, ...
        nc_info.Variables.Properties.RowNames);
    % Warn on non-existent parameters
    data_not_available = setdiff( ...
        options.parameters, ...
        nc_info.Variables.Properties.RowNames);
    for i=1:length(data_not_available)
        warning("MATLAB:cdip_request_parse_workflow", ...
                    "Data name '%s' not found.", data_not_available(i));
    end
    % Add all 2D variables, if requested
    if options.all_2D_variables == true
        data_to_query = union(data_to_query, data_2D_names);
    end
    % Add 'waveFrequency' if there's any 2D data queried
    if any(ismember(data_to_query, data_2D_names))
        data_to_query = union(data_to_query, 'waveFrequency');
    end
    % Add timestamps for each data group
    groups_in_data = data_groups(data_to_query, DATA_GROUPS);
    groups_in_data = setdiff(groups_in_data, 'meta'); % omit 'meta'
    data_to_query = union(data_to_query, strcat(groups_in_data, 'Time'));
else                            % else query all data except maybe 2D data
    data_to_query = nc_info.Variables.Name';
    if options.all_2D_variables == false
        data_to_query = setdiff(data_to_query, data_2D_names);  % remove 2D
    end
end
data_to_query = sort(data_to_query);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = ncinfo_autoretry(source, name)
%NCINFO_AUTORETRY Query info and auto retry until success or max tries met
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAX_RETRIES = 5;                         % number of query retries if error

% Submit query and get info
for i = 0:MAX_RETRIES
    try
        switch nargin
            case 1
                info = ncinfo(source);
            case 2
                info = ncinfo(source, name);
            otherwise
                MException('MATLAB:cdip_request_parse_workflow:ncinfo_autoretry ', ...
                    'Invalid number of arguments');
        end
        break;
    catch ME
        if i == MAX_RETRIES
            rethrow(ME)
        else
            pause(0.5);   % pause(seconds) and retry query
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = ncread_autoretry(source, varname, start, count, stride)
%NCREAD_AUTORETRY Query data and auto retry until success or max tries met
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAX_RETRIES = 5;                         % number of query retries if error

% Submit query and get info
for i = 0:MAX_RETRIES
    try
        switch nargin
            case 2
%                 data = ncread(source, varname);
                data = ncread_lowlevel(source, varname);
            case 4
%                 data = ncread(source, varname, start, count);
                data = ncread_lowlevel(source, varname, start, count);
            case 5
%                 data = ncread(source, varname, start, count, stride);
                data = ncread_lowlevel(source, varname, start, count, stride);
            otherwise
                MException('MATLAB:cdip_request_parse_workflow:ncread_autoretry ', ...
                    'Invalid number of arguments');
        end
        break;
    catch ME
        if i == MAX_RETRIES
            rethrow(ME)
        elseif ME.identifier == "MATLAB:imagesci:netcdf:unknownLocation" || ...
                contains(ME.message, "Variable not found")
            data = NaN;
            break;          % no need to retry
        else
            pause(0.5);     % pause(seconds) and retry query
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = ncread_lowlevel(source, varname, start, count, stride)
%NCREAD_LOWLEVEL Query data using low-level NetCDF functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncid = netcdf.open(source);

try
    varid = netcdf.inqVarID(ncid, varname);

    % Replace any inf's for reading to end with actual counts
    if nargin > 2 && any(isinf(count))
        [~, ~, dimids, ~] = netcdf.inqVar(ncid,varid);
        for i=1:length(count)
            if isinf(count(i))
                [~, dimlen] = netcdf.inqDim(ncid, dimids(i));
                count(i) = dimlen - start(i) + 1;
            end
        end
    end

    % Submit query and get info
    switch nargin
        case 2
            data = netcdf.getVar(ncid, varid);
        case 4
            data = ncread(source, varname, start, count);
        case 5
            data = ncread(source, varname, start, count, stride);
        otherwise
            MException('MATLAB:cdip_request_parse_workflow:ncread_autoretry ', ...
                'Invalid number of arguments');
    end

    netcdf.close(ncid);
catch ME
    netcdf.close(ncid);
    rethrow(ME)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function datetimes = start_end_datetimes(options)
%START_END_DATETIMES Creates structure of start and end datetimes to query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datetimes.start = {};
datetimes.end = {};
if options.years(1) > 0     % invalid/missing value = -1
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
        url_query = get_url_query(options);
        waveTime = ncread_autoretry(url_query, 'waveTime');
    end
    % Substitute in netCDF start/end dates as needed
    if options.start_date ~= ""
        datetimes.start{1} = datetime(options.start_date, ...
                                      'InputFormat', 'yyyy-MM-dd', ...
                                      'TimeZone', 'UTC');
    else
        datetimes.start{1} = datetime(waveTime(1), ...
                                      'ConvertFrom', 'posixtime', ...
                                      'TimeZone', 'UTC');
    end
    if options.end_date ~= ""
        datetimes.end{1} = datetime(options.end_date, ...
                                    'InputFormat', 'yyyy-MM-dd', ...
                                    'TimeZone', 'UTC');
        datetimes.end{1}.Hour = 23;
        datetimes.end{1}.Minute = 59;
        datetimes.end{1}.Second = 59;
    else
        datetimes.end{1} = datetime(waveTime(end), ...
                                    'ConvertFrom', 'posixtime', ...
                                    'TimeZone', 'UTC');
    end
end
end

