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

    % Initialize status structure
    status = struct('success', true, ...
                   'warnings', {{}}, ...
                   'retry_count', 0, ...
                   'failed_points', [], ...
                   'error_message', '', ...
                   'timeout_count', 0);

    % Request management configuration
    REQUEST_TIMEOUT = 60;   % seconds
    RETRY_DELAY = 5;       % seconds between retries
    MAX_RETRIES = 3;
    RATE_LIMIT_PAUSE = 2;  % seconds between requests
    BATCH_SIZE = 10;       % number of requests before cooling down
    BATCH_COOLDOWN = 30;   % seconds to pause after batch

    % Input validation
    if ~isa(data_type,'string') && ~isa(data_type,'char')
        status.success = false;
        status.error_message = 'ERROR: data_type must be a string or char';
        data = [];
        return
    end
    if ~ismember(data_type, {'3-hour', '1-hour'})
        status.success = false;
        status.error_message = 'ERROR: data_type must be either "3-hour" or "1-hour"';
        data = [];
        return
    end
    if ~isstring(parameter)
        status.success = false;
        status.error_message = 'ERROR: parameter must be a string or string array';
        data = [];
        return
    end
    if ~isa(lat_lon,'numeric') || size(lat_lon,2) ~= 2
        status.success = false;
        status.error_message = 'ERROR: lat_lon must be a numeric array with 2 columns [lat, lon]';
        data = [];
        return
    end
    if ~isa(year, 'numeric') || year < 1979 || year > 2010
        status.success = false;
        status.error_message = 'ERROR: year must be numeric and between 1979-2010';
        data = [];
        return
    end
    if ~isa(api_key,'string') && ~isa(api_key,'char')
        status.success = false;
        status.error_message = 'ERROR: api_key must be a string or char';
        data = [];
        return
    end

    % Request queue tracking
    request_queue = struct('url', {}, 'retries', {}, 'completed', {}, 'response', {}, 'point_idx', {}, 'param_name', {});
    current_batch = 0;

    try
        % determine region of lat_lon
        region = region_selection(lat_lon, data_type);

        % set domain
        if isequal(data_type,'1-hour')
            dom = ['&domain=%2Fnrel%2FUS_wave%2Fvirtual_buoy%2F' region '%2F' region '_virtual_buoy_' num2str(year) '.h5'];
        elseif isequal(data_type,'3-hour')
            dom = ['&domain=%2Fnrel%2FUS_wave%2F' region '%2F' region '_wave_' num2str(year) '.h5'];
        end

        % Get database links with timeout handling
        baseURL = ['https://developer.nrel.gov/api/hsds/?api_key=' api_key dom];
        t_start = tic;
        options = weboptions('Timeout', min(REQUEST_TIMEOUT, 30));

        try
            root = cached_webread(baseURL, options);
            if toc(t_start) > REQUEST_TIMEOUT
                status.timeout_count = status.timeout_count + 1;
                status.warnings{end+1} = 'Initial request exceeded timeout';
                error('Initial request timeout');
            end
        catch ME
            status.success = false;
            status.error_message = 'Failed to establish initial connection';
            data = [];
            return
        end

        % Get groups with timeout
        groupsURL = ['https://developer.nrel.gov/api/hsds/groups/' root.root '/links?api_key=' api_key dom];
        t_start = tic;
        try
            groups_data = cached_webread(groupsURL, options);
            if toc(t_start) > REQUEST_TIMEOUT
                status.timeout_count = status.timeout_count + 1;
                status.warnings{end+1} = 'Groups request exceeded timeout';
                error('Groups request timeout');
            end
        catch ME
            pause(RETRY_DELAY + rand()*2);
            groups_data = cached_webread(groupsURL, options);
        end
        groups = struct2table(groups_data.links);

        % Get standard parameters
        if isequal(data_type, "1-hour")
            vars = ["time_index","meta","coordinates","frequency","direction"];
        else
            vars = ["time_index","meta","coordinates"];
        end

        standard_params = struct();
        for i=1:length(vars)
            ID = groups.id(find(strcmpi(vars(i),groups.title)));
            URL = ['https://developer.nrel.gov/api/hsds/datasets/' ID{:} '/value?api_key=' api_key dom];

            % Attempt request with retries
            success = false;
            for retry=1:MAX_RETRIES
                t_start = tic;
                try
                    temp = cached_webread(URL, options);
                    if toc(t_start) <= REQUEST_TIMEOUT
                        success = true;
                        break;
                    end
                    status.timeout_count = status.timeout_count + 1;
                catch ME
                    if retry < MAX_RETRIES
                        pause(RETRY_DELAY * retry + rand()*2);
                        status.retry_count = status.retry_count + 1;
                    end
                end
            end

            if ~success
                status.warnings{end+1} = sprintf('Failed to get parameter %s after %d attempts', vars(i), MAX_RETRIES);
                continue;
            end

            if isequal(vars{i},"time_index")
                standard_params.(vars(i)) = datetime(temp.value,'InputFormat','yyyy-MM-dd HH:mm:ssXXX',...
                    'TimeZone','UTC');
            else
                standard_params.(vars(i)) = temp.value;
            end

            % Rate limiting pause with jitter
            pause(RATE_LIMIT_PAUSE * (1 + rand()));
        end

        % Find index for each location
        failed_points = [];
        for y = 1:size(lat_lon,1)
            try
                gid = lat_lon(y,:);
                radius = [abs(standard_params.coordinates(:,1)-gid(1)) abs(standard_params.coordinates(:,2)-gid(2))];
                radius(:,3) = sqrt(radius(:,1).^2 + radius(:,2).^2);
                idx(y) = find(radius(:,3)==min(radius(:,3))); % find index of closest location

                if radius(idx(y),3) > 0.5  % threshold of 0.5 degrees
                    status.warnings{end+1} = sprintf('Point [%f, %f] is %.2f degrees from nearest available point', ...
                        gid(1), gid(2), radius(idx(y),3));
                end
            catch ME
                failed_points = [failed_points; y];
                status.warnings{end+1} = sprintf('Failed to process point [%f, %f]: %s', ...
                    lat_lon(y,1), lat_lon(y,2), ME.message);
                continue
            end
        end
        status.failed_points = failed_points;

        % Create metadata struct & get parameter data
        fns_meta = ["water_depth","latitude","longitude","distance_to_shore","timezone","jurisdiction"];
        for y=1:length(idx)
            if ismember(y, failed_points)
                continue
            end

            % Get metadata
            for x=1:length(fns_meta)
                fd = standard_params.meta{idx(y)};
                meta(y).(fns_meta{x}) = fd{x};
            end

            % Get parameter data with retry and timeout
            for z = 1:length(parameter)
                paramID = groups.id(find(strcmpi(parameter(z),groups.title)));
                if isequal(parameter(z),'directional_wave_spectrum')
                    continue
                else
                    paramURL = ['https://developer.nrel.gov/api/hsds/datasets/' paramID{:} '/value?api_key=' api_key '&select=[:,' num2str(idx(y)-1) ']' dom];

                    % Attempt request with retries
                    success = false;
                    for retry=1:MAX_RETRIES
                        t_start = tic;
                        try
                            param = cached_webread(paramURL, options);
                            if toc(t_start) <= REQUEST_TIMEOUT
                                success = true;
                                break;
                            end
                            status.timeout_count = status.timeout_count + 1;
                        catch ME
                            if retry < MAX_RETRIES
                                pause(RETRY_DELAY * retry + rand()*2);
                                status.retry_count = status.retry_count + 1;
                            end
                        end
                    end

                    if success
                        param_fieldname = strrep(parameter(z), "-", "_");
                        data(y).(param_fieldname) = param.value;
                    else
                        status.warnings{end+1} = sprintf('Failed to get parameter %s for point %d', parameter(z), y);
                    end

                    % Rate limiting pause
                    pause(RATE_LIMIT_PAUSE + rand());
                end
            end

            % Batch cooldown if needed
            if mod(y, BATCH_SIZE) == 0
                status.warnings{end+1} = sprintf('Cooling down for %d seconds...', BATCH_COOLDOWN);
                pause(BATCH_COOLDOWN);
            end
        end

        % Handle directional wave spectrum if requested
        time_length = length(standard_params.time_index);
        if contains('directional_wave_spectrum',parameter)
            for y = 1:length(idx)
                if ismember(y, failed_points)
                    continue
                end

                for z = 1:time_length/486
                    i1 = num2str((z-1)*486);
                    i2 = num2str(z*486);
                    paramURL = ['https://developer.nrel.gov/api/hsds/datasets/' paramID{:} '/value?api_key=' api_key '&select=[' i1 ':' i2 ',:,:,' num2str(idx(y)-1) ']' dom];

                    success = false;
                    for retry=1:MAX_RETRIES
                        t_start = tic;
                        try
                            param = cached_webread(paramURL, options);
                            if toc(t_start) <= REQUEST_TIMEOUT
                                success = true;
                                break;
                            end
                            status.timeout_count = status.timeout_count + 1;
                        catch ME
                            if retry < MAX_RETRIES
                                pause(RETRY_DELAY * retry + rand()*2);
                                status.retry_count = status.retry_count + 1;
                            end
                        end
                    end

                    if success
                        try
                            dws = cat(1,dws,param.value);
                        catch
                            dws = param.value;
                        end
                    else
                        status.warnings{end+1} = sprintf('Failed to get spectrum data for point %d, segment %d', y, z);
                        continue;
                    end

                    % Rate limiting pause
                    pause(RATE_LIMIT_PAUSE + rand());
                end

                if exist('dws', 'var')
                    data(y).directional_wave_spectrum = dws;
                    clear dws
                end
            end
        end

        % Create final output struct
        for i=1:length(idx)
            if ismember(i, failed_points)
                continue
            end
            data(i).time = standard_params.time_index;
            if contains('directional_wave_spectrum',parameter)
                data(i).frequency = standard_params.frequency;
                data(i).direction = standard_params.direction;
            end
            data(i).metadata = meta(i);
        end

    catch ME
        status.success = false;
        status.error_message = sprintf('Fatal error: %s', ME.message);
        data = [];
        return
    end
end
