function data = request_wpto(data_type, parameter, lat_lon, year)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Returns data from the WPTO wave hindcast hosted on AWS at the specified latitude and longitude point(s), 
%     or the closest available pont(s).
%     Visit https://registry.opendata.aws/wpto-pds-us-wave/ for more information about the dataset and available 
%     locations and years. 
%     NOTE: ...................... 
%     Please see the WPTO_hindcast_example notebook for more information.  
%     
%     Parameters
%     ----------
%         data_type : string
%             Data set type of interest
%             Options: "3-hour" "1-hour"
%         parameter: string or list of strings
%             Dataset parameter to be downloaded
%             3-hour dataset options: 'directionality_coefficient', 'energy_period', 'maximum_energy_direction'
%                 'mean_absolute_period', 'mean_zero-crossing_period', 'omni-directional_wave_power', 'peak_period'
%                 'significant_wave_height', 'spectral_width', 'water_depth' 
%             1-hour dataset options: 'directionality_coefficient', 'energy_period', 'maximum_energy_direction'
%                 'mean_absolute_period', 'mean_zero-crossing_period', 'omni-directional_wave_power', 'peak_period'
%                 'significant_wave_height', 'spectral_width', 'water_depth', 'maximim_energy_direction',
%                 'mean_wave_direction', 'frequency_bin_edges'
%         lat_lon: tuple or list of tuples
%             Latitude longitude pairs at which to extract data. 
%         year : float 
%             Year to be accessed. The years 1979-2010 available.
%      
%      Returns
%      -------
%         data : struct 
%             Data indexed by datetime with columns named for parameter and cooresponding metadata index 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(data_type,'string') & ~isa(data_type,'char')
        error('ERROR: data_type must be a string or char')
    end
    if ~isa(parameter,'string')
        error('ERROR: parameter must be a string or string array')
    end
    if ~isa(lat_lon,'numeric')
        error('ERROR: lat_lon must be a double or double array')
    end
    if ~isa(year, 'numeric')
        error('ERROR: year must be numeric')
    end


    % get region
    region = region_selection(lat_lon);
    % get URL
    if isequal(data_type,"3-hour")
        url = strcat('s3://wpto-pds-us-wave/v1.0.0/',region,'/',region,'_wave_',num2str(year),'.h5');
    elseif isequal(data_type,"1-hour")
        url = strcat('s3://wpto-pds-us-wave/v1.0.0/virtual_buoy/',region,'/',region,'_virtual_buoy_',num2str(year),'.h5');
    end

    % save AWS groups into variables
    time_all = h5read(url,'/time_index/'); % read time index dataset
    meta_all = h5read(url,'/meta/'); % read metadata dataset
    coords_all = [meta_all.latitude meta_all.longitude]; % read coordinates
    % get data for each parameter
    for y = 1:length(parameter)
        data_all.(parameter(y)) = h5read(url,strcat('/',parameter(y),'/'));
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
            if isequal(f_meta{x},'jurisdiction') % format jurisdiction field
                meta(y).(f_meta{x}) = deblank(convertCharsToStrings(meta_all.jurisdiction(:,idx(y))));
            else
                meta(y).(f_meta{x}) = meta_all.(f_meta{x})(idx(y));
            end
        end
    end

    % create final output struct for data
    for i=1:length(idx)
        data(i).metadata = meta(i);
        data(i).time = time_all;
        for j=1:length(parameter)
            data(i).(parameter(j)) = data_all.(parameter(j))(idx(i),:)';
        end
    end
end
