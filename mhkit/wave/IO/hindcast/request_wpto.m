
clear; close all; clc;
%% INPUTS
tic
%
data_type = "1-hour";
year = 2010;
lat_lon = [44.624076,-124.280097; 43.489171,-125.152137]; 
param = ["significant_wave_height","energy_period"];

%% FUNCTION

region = region_selection(lat_lon);
if isequal(data_type,"3-hour")
    url = strcat('s3://wpto-pds-us-wave/v1.0.0/',region,'/',region,'_wave_',num2str(year),'.h5');
elseif isequal(data_type,"1-hour")
    url = strcat('s3://wpto-pds-us-wave/v1.0.0/virtual_buoy/',region,'/',region,'_virtual_buoy_',num2str(year),'.h5');
end

% save AWS data into variables
time_all = h5read(url,'/time_index/'); % time index
meta_all = h5read(url,'/meta/'); % metadata
coords_all = [meta_all.latitude meta_all.longitude]; % coordinates
% get data for each param
for y = 1:length(param)
    data_all.(param(y)) = h5read(url,strcat('/',param(y),'/'));
end

% find index for each location
for y = 1:size(lat_lon,1)
    gid = lat_lon(y,:);
    radius = [abs(coords_all(:,1)-gid(1)) abs(coords_all(:,2)-gid(2))];
    radius(:,3) = sqrt(radius(:,1).^2 + radius(:,2).^2);
    idx(y) = find(radius(:,3)==min(radius(:,3)));
end

% create struct for metadata
f_meta = fieldnames(meta_all);
for y=1:length(idx)
    for x=1:length(f_meta)
        if isequal(f_meta{x},'jurisdiction')
            meta(y).(f_meta{x}) = deblank(convertCharsToStrings(meta_all.jurisdiction(:,idx(y))));
        else
            meta(y).(f_meta{x}) = meta_all.(f_meta{x})(idx(y));
        end
    end
end

% create struct for data
for i=1:length(idx)
    data(i).metadata = meta(i);
    data(i).time = time_all;
    for j=1:length(param)
        data(i).(param(j)) = data_all.(param(j))(idx(i),:)';
    end
end

toc

%%
% url = strcat('s3://wpto-pds-us-wave/v1.0.0/virtual_buoy/West_Coast/');
% rfile = strcat(region,'/',region,'_virtual_buoy_',num2str(year),'.h5');
% fds = fileDatastore(url,'ReadFcn',@AWSRead);
% subds = subset(fds,1);


%% HELPERS

function region = region_selection(lat_lon)
    regions.Hawaii.lat = [15.0,15.0,27.000002,27.000002]';
    regions.Hawaii.lon = [-164.0,-151.0,-151.0,-164.0]';
    regions.West_Coast.lat = [30.0906,30.0906, 48.8641,48.8641]';
    regions.West_Coast.lon = [-130.072,-116.899,-116.899,-130.072]';
    regions.Atlantic.lat = [24.382,24.382, 44.8247,44.8247]';
    regions.Atlantic.lon = [-81.552,-65.721, -65.721,-81.552]';

    fns = fieldnames(regions);
    for i = 1:length(fns)
        in = inpolygon(lat_lon(:,1),lat_lon(:,2),regions.(fns{i}).lat,regions.(fns{i}).lon);
        if mean(in) == 1
            region = fns{i};
        end
    end
    if exist('region','var')==0
        error("WARNING: coordinates out of bounds OR not within the same region!")
    end
end

function data = AWSRead(fileName)
    data = h5read(fileName,'/meta/');
end