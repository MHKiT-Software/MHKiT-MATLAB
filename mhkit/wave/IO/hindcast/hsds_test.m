clear; close all; clc;
%% INPUTS
tic

data_type = '1-hour';
year = 2010;
parameter = ["energy_period","significant_wave_height"];
%lat_lon = [44.624076,-124.280097];
lat_lon = [44.624076,-124.280097; 43.489171,-125.152137]; 


%% REST API setup

% get API
api = getenv('HS_API_KEY');

% determine region of lat_lon
region = region_selection(lat_lon);

% set domain
if isequal(data_type,'1-hour')
    dom = ['&domain=%2Fnrel%2FUS_wave%2Fvirtual_buoy%2F' region '%2F' region '_virtual_buoy_' num2str(year) '.h5'];
elseif isequal(data_type,'3-hour')
    dom = ['&domain=%2Fnrel%2FUS_wave%2F' region '%2F' region '_wave_' num2str(year) '.h5'];
end

% get links to the databases
baseURL = ['https://developer.nrel.gov/api/hsds/?api_key=' api dom];
root = webread(baseURL);
groupsURL = ['https://developer.nrel.gov/api/hsds/groups/' root.root '/links?api_key=' api dom];
groups = webread(groupsURL);
groups = struct2table(groups.links);

% get standard parameters
vars = ["time_index","meta","coordinates"];
options = weboptions('Timeout',30);
for i=1:length(vars)
    ID = groups.id(find(strcmpi(vars(i),groups.title)));
    URL = ['https://developer.nrel.gov/api/hsds/datasets/' ID{:} '/value?api_key=' api dom];
    temp = webread(URL,options);
    standard_params.(vars(i)) = temp.value;
end

% find index for each location
for y = 1:size(lat_lon,1)
    gid = lat_lon(y,:);
    radius = [abs(standard_params.coordinates(:,1)-gid(1)) abs(standard_params.coordinates(:,2)-gid(2))];
    radius(:,3) = sqrt(radius(:,1).^2 + radius(:,2).^2);
    idx(y) = find(radius(:,3)==min(radius(:,3)));
end

% create metadata struct & get parameter data 
fns_meta = ["water_depth","latitude","longitude","distance_to_shore","timezone","jurisdiction"];
for y=1:length(idx)
    for x=1:length(fns_meta)
        fd = standard_params.meta{idx(y)};
        meta(y).(fns_meta{x}) = fd{x};
    end
    for z = 1:length(parameter)
        paramID = groups.id(find(strcmpi(parameter(z),groups.title)));
        paramURL = ['https://developer.nrel.gov/api/hsds/datasets/' paramID{:} '/value?api_key=' api '&select=[:,' num2str(idx(y)) ']' dom];
        param = webread(paramURL,options);
        data(y).(parameter(z)) = param.value; 
    end
end

% create final output struct for data
for i=1:length(idx)
    data(i).time = standard_params.time_index;
    data(i).metadata = meta(i);
end

