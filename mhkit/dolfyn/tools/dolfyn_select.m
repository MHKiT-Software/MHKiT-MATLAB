function [selected_data, selected_coords] = dolfyn_select(data_field, coord_name, target_value, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select ADCP data along coordinates/dimensions with xarray-like functionality
%
% This function provides data selection capabilities similar to xarray's .sel() method,
% enabling efficient selection of  velocity, range, time, or beam data along
% specified coordinate dimensions.
%
% Parameters
% ------------
%   data_field : structure
%       ADCP data structure containing:
%         .data : Velocity or other data array [varies]
%         .coords : Coordinate structure with dimension arrays
%   coord_name : string
%       Name of coordinate dimension to select along
%   target_value : double
%       Value(s) to select along coordinate (scalar or array)
%   method : string, optional (name-value)
%       Selection method: 'exact', 'nearest', 'pad', 'ffill', 'backfill', 'bfill'
%       Default: 'exact'
%   tolerance : double, optional (name-value)
%       Maximum distance for nearest method [same units as coordinate]
%       Default: Inf
%   return_index : logical, optional (name-value)
%       Return indices instead of data (default: false)
%
% Returns
% ---------
%   selected_data : array or double
%       Selected data array or coordinate indices [same units as input data]
%   selected_coords : structure
%       Updated coordinate structure with selected coordinate values
%       .{coord_name} : Selected coordinate array [same units as input]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data_field (1,1) struct
    coord_name {mustBeTextScalar}
    target_value (:,1) {mustBeNumeric}
    options.method {mustBeTextScalar} = 'exact'
    options.tolerance (1,1) {mustBeNumeric, mustBePositive} = Inf
    options.return_index (1,1) logical = false
end

% Extract parameters
method = options.method;
tolerance = options.tolerance;
return_index = options.return_index;

% Validate data structure
if ~isfield(data_field, 'data')
    error('MHKiT:dolfyn_select:InvalidInput', 'data_field must contain a .data field');
end
if ~isfield(data_field, 'coords')
    error('MHKiT:dolfyn_select:InvalidInput', 'data_field must contain a .coords field');
end
if ~isfield(data_field.coords, coord_name)
    error('MHKiT:dolfyn_select:InvalidCoord', 'Coordinate "%s" not found in data_field.coords', coord_name);
end

% Get coordinate values
coord_values = data_field.coords.(coord_name);
if isempty(coord_values)
    error('MHKiT:dolfyn_select:EmptyCoord', 'Coordinate "%s" is empty', coord_name);
end

% Find coordinate dimension in data
coord_dims = fieldnames(data_field.coords);
coord_idx = find(strcmp(coord_dims, coord_name));
if isempty(coord_idx)
    error('MHKiT:dolfyn_select:CoordNotFound', 'Could not determine dimension for coordinate "%s"', coord_name);
end

% Handle different selection methods
selected_indices = [];
for i = 1:length(target_value)
    target = target_value(i);
    
    switch lower(method)
        case 'exact'
            % Exact match
            idx = find(coord_values == target);
            if isempty(idx)
                error('MHKiT:dolfyn_select:ExactMatchNotFound', ...
                    'Exact match not found for value %.6f in coordinate "%s"', target, coord_name);
            end
            
        case 'nearest'
            % Nearest value
            [min_diff, idx] = min(abs(coord_values - target));
            if min_diff > tolerance
                error('MHKiT:dolfyn_select:ToleranceExceeded', ...
                    'Nearest value %.6f exceeds tolerance %.6f for target %.6f', ...
                    min_diff, tolerance, target);
            end
            
        case {'pad', 'ffill'}
            % Forward fill - find nearest lower value
            valid_indices = find(coord_values <= target);
            if isempty(valid_indices)
                error('MHKiT:dolfyn_select:NoValidPadValue', ...
                    'No coordinate values <= %.6f found for pad/ffill', target);
            end
            [~, max_idx] = max(coord_values(valid_indices));
            idx = valid_indices(max_idx);
            
        case {'backfill', 'bfill'}
            % Backward fill - find nearest higher value  
            valid_indices = find(coord_values >= target);
            if isempty(valid_indices)
                error('MHKiT:dolfyn_select:NoValidBackfillValue', ...
                    'No coordinate values >= %.6f found for backfill/bfill', target);
            end
            [~, min_idx] = min(coord_values(valid_indices));
            idx = valid_indices(min_idx);
            
        otherwise
            error('MHKiT:dolfyn_select:InvalidMethod', ...
                'Invalid method "%s". Valid methods: exact, nearest, pad, ffill, backfill, bfill', method);
    end
    
    selected_indices = [selected_indices, idx];
end

% Handle multiple indices (remove duplicates, preserve order)
[selected_indices, ~] = unique(selected_indices, 'stable');

% Return indices if requested
if return_index
    selected_data = selected_indices;
    selected_coords = data_field.coords;
    return;
end

% Select data based on coordinate dimension
data_size = size(data_field.data);
ndims_data = length(data_size);

% Create indexing cell array
indices = cell(1, ndims_data);
for dim = 1:ndims_data
    indices{dim} = ':';
end

% Determine which dimension corresponds to the coordinate
% This is a simplified heuristic - in practice, you might need more sophisticated logic
if length(coord_values) == data_size(coord_idx)
    % Coordinate dimension matches data dimension
    indices{coord_idx} = selected_indices;
elseif coord_idx <= ndims_data && length(coord_values) == data_size(coord_idx)
    indices{coord_idx} = selected_indices;
else
    % Try to find matching dimension
    matching_dim = find(data_size == length(coord_values));
    if length(matching_dim) == 1
        indices{matching_dim} = selected_indices;
    else
        % Default to last dimension (common in dolfyn)
        indices{end} = selected_indices;
    end
end

% Extract selected data
selected_data = data_field.data(indices{:});

% Update coordinates
selected_coords = data_field.coords;
selected_coords.(coord_name) = coord_values(selected_indices);

end
