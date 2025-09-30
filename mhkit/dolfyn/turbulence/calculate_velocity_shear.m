function ds_out = calculate_velocity_shear(ds, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate velocity shear components (dudz, dvdz, dwdz) from ADCP velocity profiles.
%
% Parameters
% ------------
%   ds: structure
%       ADCP dataset structure containing velocity data
%   vel_field: string
%       Name of the velocity field in the dataset
%       Default: 'vel'
%   diff_style: string
%       Differentiation method: 'first', 'centered', or 'centered_extended'
%       Default: 'centered_extended'
%   components: cell array
%       Velocity components to calculate shear for: {'u', 'v', 'w'}
%       Default: {'u', 'v', 'w'} (all components)
%   calc_shear_squared: logical
%       Calculate horizontal shear squared (dudz² + dvdz²)
%       Default: true
%   field_names: cell array
%       Names of output fields: {dudz, dvdz, dwdz, shear_squared}
%       Default: {'dudz', 'dvdz', 'dwdz', 'shear_squared'}
%
% Returns
% ---------
%   ds_out: structure
%       Input dataset with added shear fields:
%           ds_out.dudz.data : du/dz shear values [s⁻¹]
%           ds_out.dvdz.data : dv/dz shear values [s⁻¹] 
%           ds_out.dwdz.data : dw/dz shear values [s⁻¹]
%           ds_out.shear_squared.data : Horizontal shear squared [s⁻²]
%           All fields include proper dims, coords, units, and metadata
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds
        options.vel_field = 'vel'
        options.diff_style = 'centered_extended'
        options.components = {'u', 'v', 'w'}
        options.calc_shear_squared = true
        options.field_names = {'dudz', 'dvdz', 'dwdz', 'shear_squared'}
    end
    
    % Validate input dataset structure
    if ~isstruct(ds)
        error('mhkit:dolfyn:calculate_velocity_shear: Input ds must be a structure');
    end
    
    % Check for required velocity field
    if ~isfield(ds, options.vel_field)
        error('mhkit:dolfyn:calculate_velocity_shear: Dataset must contain velocity field: %s', options.vel_field);
    end
    
    vel_data = ds.(options.vel_field);
    
    % Validate velocity structure
    if ~isfield(vel_data, 'data') || ~isfield(vel_data, 'dims') || ~isfield(vel_data, 'coords')
        error('mhkit:dolfyn:calculate_velocity_shear: Velocity field must contain data, dims, and coords');
    end
    
    vel_values = vel_data.data;
    
    % Handle data with singleton dimensions by squeezing them out
    original_size = size(vel_values);
    vel_values = squeeze(vel_values);
    
    % Validate velocity data dimensions after squeezing
    if ndims(vel_values) ~= 3
        error('mhkit:dolfyn:calculate_velocity_shear: Velocity data must be 3D (range x time x dir) after removing singleton dimensions. Got size: %s', ...
            mat2str(original_size));
    end
    
    % Find velocity component and range dimensions
    dir_dim = find(strcmp(vel_data.dims, 'dir'));
    range_dim = find(strcmp(vel_data.dims, 'range'));
    time_dim = find(strcmp(vel_data.dims, 'time'));
    
    if isempty(dir_dim) || isempty(range_dim) || isempty(time_dim)
        error('mhkit:dolfyn:calculate_velocity_shear: Velocity data must have range, time, and dir dimensions');
    end
    
    % Ensure velocity data is in [range, time, dir] order
    if range_dim ~= 1 || time_dim ~= 2 || dir_dim ~= 3
        % Reorder dimensions to [range, time, dir]
        perm_order = [range_dim, time_dim, dir_dim];
        vel_values = permute(vel_values, perm_order);
    end
    
    [n_range, n_time, n_dir] = size(vel_values);
    
    % Check minimum range requirement
    if n_range < 2
        error('mhkit:dolfyn:calculate_velocity_shear: Need at least 2 range bins to calculate shear');
    end
    
    % Get range coordinate
    if isfield(vel_data.coords, 'range')
        range_coord = vel_data.coords.range;
        range_name = 'range';
    else
        error('mhkit:dolfyn:calculate_velocity_shear: Velocity data must contain range coordinate');
    end
    
    % Get time coordinate
    if isfield(vel_data.coords, 'time')
        time_coord = vel_data.coords.time;
        time_name = 'time';
    else
        error('mhkit:dolfyn:calculate_velocity_shear: Velocity data must contain time coordinate');
    end
    
    % Check minimum number of velocity components
    if n_dir < 3 && any(strcmp(options.components, 'w'))
        warning('mhkit:dolfyn:calculate_velocity_shear: Dataset has only %d velocity components, skipping w-component', n_dir);
        options.components = options.components(~strcmp(options.components, 'w'));
    end
    
    % Validate differentiation style
    valid_styles = {'first', 'centered', 'centered_extended'};
    if ~ismember(options.diff_style, valid_styles)
        error('mhkit:dolfyn:calculate_velocity_shear: diff_style must be one of: %s', strjoin(valid_styles, ', '));
    end

    % Initialize output
    ds_out = ds;
    
    % Component mapping
    comp_map = containers.Map({'u', 'v', 'w'}, {1, 2, 3});
    field_map = containers.Map({'u', 'v', 'w'}, options.field_names(1:3));
    
    % Calculate shear for each requested component
    for i = 1:length(options.components)
        comp = options.components{i};
        
        if ~isKey(comp_map, comp)
            warning('Unknown velocity component: %s', comp);
            continue;
        end
        
        comp_idx = comp_map(comp);
        if comp_idx > n_dir
            warning('Component %s not available (only %d components in data)', comp, n_dir);
            continue;
        end
        
        field_name = field_map(comp);
        
        % Extract velocity component
        vel_comp = squeeze(vel_values(:, :, comp_idx));  % [range x time]
        
        % Calculate shear using specified method
        [shear_values, range_out] = calculate_shear_profile(vel_comp, range_coord, options.diff_style);
        
        % Create output field
        ds_out.(field_name) = struct();
        ds_out.(field_name).data = single(shear_values);
        ds_out.(field_name).dims = {range_name, time_name};
        ds_out.(field_name).coords = struct();
        ds_out.(field_name).coords.(range_name) = range_out;
        ds_out.(field_name).coords.(time_name) = time_coord;
        ds_out.(field_name).units = "s^-1";
        ds_out.(field_name).long_name = sprintf("Shear in %s-direction", upper(comp));
        ds_out.(field_name).method = options.diff_style;
        ds_out.(field_name).description = sprintf('d%s/dz calculated using %s difference', comp, options.diff_style);
    end
    
    % Calculate horizontal shear squared (dudz² + dvdz²)
    if options.calc_shear_squared && isfield(ds_out, options.field_names{1}) && isfield(ds_out, options.field_names{2})
        
        dudz_data = ds_out.(options.field_names{1}).data;
        dvdz_data = ds_out.(options.field_names{2}).data;
        
        shear_squared = dudz_data.^2 + dvdz_data.^2;
        
        % Create shear squared field
        ds_out.(options.field_names{4}) = struct();
        ds_out.(options.field_names{4}).data = single(shear_squared);
        ds_out.(options.field_names{4}).dims = {range_name, time_name};
        ds_out.(options.field_names{4}).coords = ds_out.(options.field_names{1}).coords;  % Same coords as dudz
        ds_out.(options.field_names{4}).units = "s^-2";
        ds_out.(options.field_names{4}).long_name = "Horizontal Shear Squared";
        ds_out.(options.field_names{4}).description = sprintf('(du/dz)² + (dv/dz)² calculated using %s method', options.diff_style);
    end

end

function [shear_out, range_out] = calculate_shear_profile(vel_profile, range_coord, diff_style)
    % Calculate velocity shear using specified differentiation method
    
    n_range = size(vel_profile, 1);
    n_time = size(vel_profile, 2);
    
    % Calculate range differences
    dz = diff(range_coord);
    
    switch diff_style
        case 'first'
            % Forward difference: (u[i+1] - u[i]) / (z[i+1] - z[i])
            dz_column = dz(:);  % Ensure column vector for proper broadcasting
            shear_out = diff(vel_profile, 1, 1) ./ dz_column;  % [n_range-1 x n_time]
            range_out = range_coord(2:end);  % Drop first range point
            
        case 'centered'
            % Centered difference: (u[i+1] - u[i-1]) / (2 * dz)
            if n_range < 3
                error('mhkit:dolfyn:calculate_shear_profile: Centered difference requires at least 3 range bins');
            end
            
            % Calculate centered difference
            vel_forward = vel_profile(3:end, :);      % u[i+1]  [n_range-2 x n_time]
            vel_backward = vel_profile(1:end-2, :);   % u[i-1]  [n_range-2 x n_time]
            dz_centered = 2 * dz(2:end);             % 2 * (z[i+1] - z[i])  [n_range-2 x 1]
            
            % Ensure proper broadcasting by making dz_centered a column vector
            dz_centered = dz_centered(:);  % [n_range-2 x 1]
            shear_out = (vel_forward - vel_backward) ./ dz_centered;  % [n_range-2 x n_time]
            range_out = range_coord(2:end-1);  % Drop first and last range points
            
        case 'centered_extended'
            % Centered difference with edge preservation
            if n_range < 3
                % Fall back to first difference if insufficient points
                [shear_out, range_out] = calculate_shear_profile(vel_profile, range_coord, 'first');
                return;
            end
            
            % Initialize output
            shear_out = zeros(n_range, n_time);
            
            % First point: forward difference
            shear_out(1, :) = (vel_profile(2, :) - vel_profile(1, :)) / dz(1);
            
            % Middle points: centered difference
            for i = 2:n_range-1
                shear_out(i, :) = (vel_profile(i+1, :) - vel_profile(i-1, :)) / (2 * dz(i));
            end
            
            % Last point: backward difference  
            shear_out(end, :) = (vel_profile(end, :) - vel_profile(end-1, :)) / dz(end);
            
            range_out = range_coord;  % Keep all range points
            
        otherwise
            error('mhkit:dolfyn:calculate_shear_profile: Unknown differentiation style: %s', diff_style);
    end
end
