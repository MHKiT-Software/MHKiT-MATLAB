function data = request_wpto(data_type, parameter, lat_lon, year, api_key)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Returns data from the WPTO wave hindcast hosted on AWS at the specified latitude and longitude point(s), 
%     or the closest available pont(s).
%     Visit https://registry.opendata.aws/wpto-pds-us-wave/ for more information about the dataset and available 
%     locations and years. 
%     NOTE: To access the WPTO hindcast data, you will need to configure h5pyd for data access on HSDS. 
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
%                 'mean_wave_direction', 'frequency_bin_edges', 'directional_wave_spectrum'
%         lat_lon: tuple or list of tuples
%             Latitude longitude pairs at which to extract data. 
%         year : float 
%             Year to be accessed. The years 1979-2010 available.
%         api_key : string
%             API key obtained from https://developer.nrel.gov/signup/
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
    if ~isstring(parameter)
        error('ERROR: parameter must be a string or string array')
    end
    if ~isa(lat_lon,'numeric')
        error('ERROR: lat_lon must be a double or double array')
    end
    if ~isa(year, 'numeric')
        error('ERROR: year must be numeric')
    end
    if ~isa(api_key,'string') & ~isa(api_key,'char')
        error('ERROR: api_key must be a string or string array')
    end


    % set weboptions
    options = weboptions('Timeout',30);

    % determine region of lat_lon
    region = region_selection(lat_lon);

    % set domain
    if isequal(data_type,'1-hour')
        dom = ['&domain=%2Fnrel%2FUS_wave%2Fvirtual_buoy%2F' region '%2F' region '_virtual_buoy_' num2str(year) '.h5'];
    elseif isequal(data_type,'3-hour')
        dom = ['&domain=%2Fnrel%2FUS_wave%2F' region '%2F' region '_wave_' num2str(year) '.h5'];
    end

    % get links to the databases
    baseURL = ['https://developer.nrel.gov/api/hsds/?api_key=' api_key dom];
    root = webread(baseURL);
    groupsURL = ['https://developer.nrel.gov/api/hsds/groups/' root.root '/links?api_key=' api_key dom];
    groups = webread(groupsURL);
    groups = struct2table(groups.links);

    % get standard parameters
    if isequal(data_type, "1-hour")
        vars = ["time_index","meta","coordinates","frequency","direction"];
    else
        vars = ["time_index","meta","coordinates"];
    end
    for i=1:length(vars)
        ID = groups.id(find(strcmpi(vars(i),groups.title)));
        URL = ['https://developer.nrel.gov/api/hsds/datasets/' ID{:} '/value?api_key=' api_key dom];
        temp = webread(URL,options);
        standard_params.(vars(i)) = temp.value;
    end

    % find index for each location
    for y = 1:size(lat_lon,1)
        gid = lat_lon(y,:);
        radius = [abs(standard_params.coordinates(:,1)-gid(1)) abs(standard_params.coordinates(:,2)-gid(2))];
        radius(:,3) = sqrt(radius(:,1).^2 + radius(:,2).^2);
        idx(y) = find(radius(:,3)==min(radius(:,3))); % find index of closest location
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
            if isequal(parameter(z),'directional_wave_spectrum')
                continue
            else
                paramURL = ['https://developer.nrel.gov/api/hsds/datasets/' paramID{:} '/value?api_key=' api_key '&select=[:,' num2str(idx(y)-1) ']' dom];
            end
            param = webread(paramURL,options);
            data(y).(parameter(z)) = param.value; 
        end
    end

    % handle extraction of directional wave spectrum if requested
    time_length = length(standard_params.time_index);
    if contains('directional_wave_spectrum',parameter)
        for y = 1:length(idx)
            for z = 1:time_length/486
                i1 = num2str((z-1)*486);
                i2 = num2str(z*486);
                paramURL = ['https://developer.nrel.gov/api/hsds/datasets/' paramID{:} '/value?api_key=' api_key '&select=[' i1 ':' i2 ',:,:,' num2str(idx(y)-1) ']' dom];
                param = webread(paramURL,options);
                try
                    dws = cat(1,dws,param.value); 
                catch
                    dws = param.value;
                end
            end
            data(y).directional_wave_spectrum = dws;
            clear dws
        end  
    end

    % create final output struct for data
    for i=1:length(idx)
        data(i).time = standard_params.time_index;
        if contains('directional_wave_spectrum',parameter)
            data(i).frequency = standard_params.frequency;
            data(i).direction = standard_params.direction;
        end
        data(i).metadata = meta(i);
    end
end