function [data] = cdip_request_parse_workflow(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Request and parse CDIP buoy data from the CDIP THREDDS server
%
% Requests data for a station number (from http://cdip.ucsd.edu/) and
% parses it into grouped structures. Years may be non-consecutive, e.g.
% [2001, 2010]. Time may be sliced by dates (start_date or end_date in
% YYYY-MM-DD). 2D variables are only returned when requested.
%
% Parameters
% ------------
%     station_number : string
%         Station number of the CDIP wave buoy
%     parameters : string or array of strings (optional)
%         Variables to return. Default returns all variables except 2D
%         variables
%     years : int or array of int (optional)
%         Year, e.g. 2001 or [2001, 2010]
%     start_date : string (optional)
%         Start date in YYYY-MM-DD, e.g. '2012-04-01'
%     end_date : string (optional)
%         End date in YYYY-MM-DD, e.g. '2012-04-30'
%     data_type : string (optional)
%         'historic' (default) or 'realtime'
%     all_2D_variables : logical (optional)
%         Return all 2D data. This adds significant processing time, so
%         pass the 2D variables of interest in parameters instead where
%         possible. Default false
%
% Returns
% ---------
%     data : structure
%         data.data : structure
%             Grouped 1D and 2D structures of array data with datetimes
%         data.metadata : structure
%             Anything not of length time, including the buoy name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Open the remote dataset once and read every variable through this handle.
% Opening per variable costs several HTTP requests each and gets the client
% rate limited by the CDIP THREDDS server.
ncid = netcdf_open_autoretry(url_query);
close_dataset = onCleanup(@() netcdf.close(ncid));

% Build list of data to query
data_to_query = make_data_list(options, nc_info, DATA_GROUPS);

% Create list of start and end datetimes/indices for which to query data
datetimes = start_end_datetimes(options, ncid);
indices = data_indices(ncid, datetimes, data_to_query, DATA_GROUPS);

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
            value = ncread_autoretry(ncid, name);
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
                    value = ncread_autoretry(ncid, name, ...
                                   [1, index_start], [Inf, index_count]);
                    value = value';
                elseif shape == "1D"
                    value = ncread_autoretry(ncid, name, ...
                                   index_start, index_count);
                end
            catch ME
                if ME.identifier == "MHKiT:cdip_request_parse_workflow:AccessDenied"
                    rethrow(ME)
                elseif ME.identifier == "MATLAB:imagesci:netcdf:libraryFailure"
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
    ncread_autoretry(ncid, 'metaStationName')));
end


function group = data_group(data_name, all_groups)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return the data group of a variable from its name prefix
%
% Parameters
% ------------
%     data_name : char
%         Variable name
%     all_groups : cell
%         Group prefixes to match
%
% Returns
% ---------
%     group : char
%         Matching group, or 'other' if no prefix matches
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

group = 'other';
for i = 1:length(all_groups)
    if startsWith(data_name, all_groups{i})     % group is the prefix
        group = all_groups{i};
        break
    end
end
end


function groups = data_groups(data_names, all_groups)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return the unique set of data groups for a list of variable names
%
% Parameters
% ------------
%     data_names : cell
%         Variable names
%     all_groups : cell
%         Group prefixes to match
%
% Returns
% ---------
%     groups : cell
%         Unique groups present in data_names
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrapper_fun = @(x) data_group(x, all_groups);
group_of_each = cellfun(wrapper_fun, data_names, ...
                           'UniformOutput', false);
groups = unique(group_of_each);
end


function indices = data_indices(ncid, datetime_ranges, data_to_query, all_groups)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return the start and end indices to query for each data group and
% time range
%
% Parameters
% ------------
%     ncid : double
%         Handle from netcdf.open
%     datetime_ranges : structure
%         datetime_ranges.start and .end cell arrays of datetimes
%     data_to_query : cell
%         Variable names to be queried
%     all_groups : cell
%         Group prefixes to match
%
% Returns
% ---------
%     indices : structure
%         indices.<group>.start and indices.<group>.end index vectors, one
%         entry per time range that holds data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groups_in_data = data_groups(data_to_query, all_groups);

indices = struct;
for i = 1:length(groups_in_data)
    posixtimes = ncread_autoretry( ...
            ncid, strcat(groups_in_data{i}, 'Time'));
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


function [type, group, shape] = datum_categories(datum_name, nc_info, all_groups)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return the type, group, and shape of a variable
%
% Parameters
% ------------
%     datum_name : char
%         Variable name
%     nc_info : structure
%         ncinfo output with Variables as a table
%     all_groups : cell
%         Group prefixes to match
%
% Returns
% ---------
%     type : char
%         'data' if the variable is of length time, else 'metadata'
%     group : char
%         'wave', 'sst', 'gps', 'dwr', 'meta', or 'other'
%     shape : char
%         '0D', '1D', or '2D'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function data_2D_names = find_data_2D_names(nc_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find the names of the 2D (frequency by time) variables
%
% Parameters
% ------------
%     nc_info : structure
%         ncinfo output with Variables as a table
%
% Returns
% ---------
%     data_2D_names : cell
%         Variable names sized [waveFrequency, waveTime]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_freq = nc_info.Variables{'waveFrequency', 'Size'}{1};
N_time = nc_info.Variables{'waveTime', 'Size'}{1};
data_2D_names = {};
for i = 1:height(nc_info.Variables)
    if isequal(nc_info.Variables.Size{i}, [N_freq, N_time])
        data_2D_names{end+1} = nc_info.Variables.Name{i};
    end
end
end


function url_query = get_url_query(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Build the OPeNDAP URL of the station dataset
%
% Parameters
% ------------
%     options : structure
%         Parsed name-value options with station_number and data_type
%
% Returns
% ---------
%     url_query : string
%         THREDDS URL of the historic or realtime NetCDF file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function data_to_query = make_data_list(options, nc_info, DATA_GROUPS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compile the sorted list of variables to query
%
% Requested parameters are checked against the dataset, and the
% waveFrequency and <group>Time variables they depend on are added.
%
% Parameters
% ------------
%     options : structure
%         Parsed name-value options
%     nc_info : structure
%         ncinfo output with Variables as a table
%     DATA_GROUPS : cell
%         Group prefixes
%
% Returns
% ---------
%     data_to_query : cell
%         Sorted variable names
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function info = ncinfo_autoretry(source)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Query dataset info, retrying transient failures
%
% Parameters
% ------------
%     source : string
%         Dataset URL
%
% Returns
% ---------
%     info : structure
%         ncinfo output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = retry_remote_request(@() ncinfo(source));
end


function ncid = netcdf_open_autoretry(source)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Open the remote dataset read only, retrying transient failures
%
% Parameters
% ------------
%     source : string
%         Dataset URL
%
% Returns
% ---------
%     ncid : double
%         Handle from netcdf.open
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncid = retry_remote_request(@() netcdf.open(source, 'NOWRITE'));
end


function data = ncread_autoretry(ncid, varname, start, count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read a variable, or a slice of it, retrying transient failures
%
% Parameters
% ------------
%     ncid : double
%         Handle from netcdf.open
%     varname : char
%         Variable name
%     start : vector (optional)
%         One-based start index per dimension
%     count : vector (optional)
%         Elements to read per dimension, Inf reads to the end
%
% Returns
% ---------
%     data : array
%         Variable data, or NaN if the variable does not exist
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    request = @() read_variable(ncid, varname);
else
    request = @() read_variable(ncid, varname, start, count);
end
try
    data = retry_remote_request(request);
catch ME
    if contains(ME.message, "Variable not found")
        data = NaN;
    else
        rethrow(ME)
    end
end
end


function result = retry_remote_request(request)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run a remote NetCDF request, retrying transient failures with
% exponential backoff
%
% Missing variables and an access denial from the server are not
% transient, so they are raised immediately. Retrying after a denial
% only extends the block the CDIP THREDDS server has placed on the client.
%
% Parameters
% ------------
%     request : function handle
%         Zero-argument function that performs the request
%
% Returns
% ---------
%     result : any
%         Whatever the request returns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAX_RETRIES = 5;
for i = 0:MAX_RETRIES
    try
        result = request();
        return
    catch ME
        if contains(ME.message, "Authorization failure")
            error('MHKiT:cdip_request_parse_workflow:AccessDenied', ...
                ['The CDIP THREDDS server refused the request (HTTP access denied). ' ...
                 'This is usually rate limiting; wait before retrying. ' ...
                 'NetCDF reported: %s'], ME.message);
        elseif contains(ME.message, "Variable not found") || i == MAX_RETRIES
            rethrow(ME)
        end
        pause(0.5 * 2^i);   % 0.5, 1, 2, 4, 8 seconds
    end
end
end


function data = read_variable(ncid, varname, start, count)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read a variable, or a slice of it, through an open NetCDF handle
%
% The sliced read applies the _FillValue, scale_factor, and add_offset
% attributes the way ncread does, so the output matches the previous
% ncread based implementation.
%
% Parameters
% ------------
%     ncid : double
%         Handle from netcdf.open
%     varname : char
%         Variable name
%     start : vector (optional)
%         One-based start index per dimension, as for ncread
%     count : vector (optional)
%         Elements to read per dimension, Inf reads to the end
%
% Returns
% ---------
%     data : array
%         Variable data in MATLAB dimension order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varid = netcdf.inqVarID(ncid, varname);
if nargin == 2
    data = netcdf.getVar(ncid, varid);
    return
end

% Replace any inf's for reading to end with actual counts
[~, ~, dimids] = netcdf.inqVar(ncid, varid);
for i = find(isinf(count))
    [~, dimlen] = netcdf.inqDim(ncid, dimids(i));
    count(i) = dimlen - start(i) + 1;
end

% netcdf.getVar takes zero-based start indices
data = netcdf.getVar(ncid, varid, start - 1, count);
data = apply_cf_attributes(ncid, varid, data);
end


function data = apply_cf_attributes(ncid, varid, data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply the _FillValue, scale_factor, and add_offset attributes to raw
% variable data following ncread
%
% Fill values become NaN in floating point output. scale_factor and
% add_offset convert the data to double.
%
% Parameters
% ------------
%     ncid : double
%         Handle from netcdf.open
%     varid : double
%         Variable id
%     data : array
%         Raw data from netcdf.getVar
%
% Returns
% ---------
%     data : array
%         Converted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fill_value = get_attribute(ncid, varid, '_FillValue');
scale_factor = get_attribute(ncid, varid, 'scale_factor');
add_offset = get_attribute(ncid, varid, 'add_offset');

is_fill = false(size(data));
if ~isempty(fill_value)
    is_fill = data == fill_value;
end
if ~isempty(scale_factor) || ~isempty(add_offset)
    data = double(data);
    if ~isempty(scale_factor)
        data = data * double(scale_factor);
    end
    if ~isempty(add_offset)
        data = data + double(add_offset);
    end
end
if isfloat(data)
    data(is_fill) = NaN;
end
end


function value = get_attribute(ncid, varid, name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return a variable attribute, or empty if it is not defined
%
% Parameters
% ------------
%     ncid : double
%         Handle from netcdf.open
%     varid : double
%         Variable id
%     name : char
%         Attribute name
%
% Returns
% ---------
%     value : any
%         Attribute value, or [] when the attribute is absent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    value = netcdf.getAtt(ncid, varid, name);
catch
    value = [];
end
end


function datetimes = start_end_datetimes(options, ncid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the start and end datetimes of each time range to query
%
% Years give one range per year. Otherwise a single range is built from
% start_date and end_date, falling back to the first and last waveTime
% in the dataset for whichever is not given.
%
% Parameters
% ------------
%     options : structure
%         Parsed name-value options
%     ncid : double
%         Handle from netcdf.open
%
% Returns
% ---------
%     datetimes : structure
%         datetimes.start and datetimes.end cell arrays of datetimes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if options.start_date == "" || options.end_date == ""
        waveTime = ncread_autoretry(ncid, 'waveTime');
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

