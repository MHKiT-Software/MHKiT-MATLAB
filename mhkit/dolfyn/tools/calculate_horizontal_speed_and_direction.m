function ds_out = calculate_horizontal_speed_and_direction(ds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate horizontal speed and direction from velocity components
%
% Parameters
% ------------
%     ds: structure
%         Dolfyn data structure containing velocity data. Must include:
%             - vel: structure with fields:
%                 - data: array with velocity components
%                 - dims: cell array of dimension names
%                 - coords: structure with coordinate information including:
%                     - dir: cell array of direction labels ('E', 'N', 'U1', 'U2')
%             - attrs: structure with metadata
%             - coord_sys: string
%                 Coordinate system identifier (e.g., "earth", "beam", "inst", "principal")
%
% Returns
% ------------
%     ds_out: structure
%         Input structure with additional fields:
%             - U_mag: structure
%                 - data: array of velocity magnitudes (m/s)
%                 - dims: cell array of dimension names (excluding dir)
%                 - coords: coordinate information
%                 - units: "m s-1"
%                 - standard_name: "sea_water_speed"
%                 - long_name: "Water Speed"
%             - U_dir: structure
%                 - data: array of flow directions
%                 - dims: cell array of dimension names (excluding dir)
%                 - coords: coordinate information
%                 - units: "degrees_CW_from_[REF]" where [REF] depends on coord_sys:
%                     - earth: "degrees_CW_from_N"
%                     - others: "degrees_CW_from_E" (or relevant streamwise component)
%                 - standard_name: "sea_water_to_direction"
%                 - long_name: "Water Direction"
%                 Note: Direction is NaN where velocity magnitude is zero
%                 Note: For earth coordinates, angles are CW from North [0, 360]
%                       For other coordinates, angles are CW from East/streamwise [0, 360]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input validation - check for required fields
    validateattributes(ds, {'struct'}, {'nonempty'}, mfilename, 'ds');
    assert(isfield(ds, 'attrs'), 'MATLAB:assert:failed', 'Field attrs missing from ds');
    assert(isfield(ds, 'vel'), 'MATLAB:assert:failed', 'Field vel missing from ds');
    validateattributes(ds.vel, {'struct'}, {'nonempty'}, mfilename, 'ds.vel');

    % Check for required vel subfields
    required_fields = {'data', 'dims', 'coords'};
    for i = 1:length(required_fields)
        assert(isfield(ds.vel, required_fields{i}), 'MATLAB:assert:failed', ...
            ['Field ' required_fields{i} ' missing from ds.vel']);
    end

    % Check for coords.dir
    assert(isfield(ds.vel.coords, 'dir'), 'MATLAB:assert:failed', ...
        'Field dir missing from ds.vel.coords');
    validateattributes(ds.vel.coords.dir, {'cell'}, {'nonempty'}, mfilename, 'ds.vel.coords.dir');

    % Find the index of the 'dir' dimension
    dir_idx = find(strcmp(ds.vel.dims, 'dir'));
    assert(~isempty(dir_idx), 'MATLAB:assert:failed', ...
        "'dir' dimension not found in vel.dims");

    % Get the size of the data array
    data_size = size(ds.vel.data);

    % Ensure data array has enough dimensions
    assert(length(data_size) >= dir_idx, 'MATLAB:assert:failed', ...
        'Data array has fewer dimensions than the dir index');

    % Find indices for East and North components in the direction coordinates
    dir_coords = ds.vel.coords.dir;

    E_idx = find(strcmp(dir_coords, 'E'));
    N_idx = find(strcmp(dir_coords, 'N'));

    % Verify both components exist
    assert(~isempty(E_idx), 'MATLAB:assert:failed', ...
        'East component not found in vel.coords.dir');
    assert(~isempty(N_idx), 'MATLAB:assert:failed', ...
        'North component not found in vel.coords.dir');

    % Create indexing cell array for the entire data array
    idx_cell = repmat({':'}, 1, length(data_size));

    % Extract East and North velocities using dynamic indexing
    idx_cell{dir_idx} = E_idx;

    % Remove singleton dimensions
    vel = squeeze(ds.vel.data);
    u = squeeze(vel(idx_cell{:}));  % East component

    idx_cell{dir_idx} = N_idx;
    v = squeeze(vel(idx_cell{:}));  % North component

    % Calculate complex velocity
    U = u + v * 1j;
    U_mag = abs(U);

    % Prepare output structure for magnitude
    ds_out = ds;
    ds_out.U_mag = struct();
    ds_out.U_mag.data = U_mag;
    ds_out.U_mag.dims = ds.vel.dims(1:end ~= dir_idx);
    ds_out.U_mag.coords = ds.vel.coords;
    ds_out.U_mag.units = "m s-1";
    ds_out.U_mag.standard_name = "sea_water_speed";
    ds_out.U_mag.long_name = "Water Speed";

    % Calculate direction
    angle_rad = angle(U);
    angle_deg = angle_rad * (180/pi);

    % Convert angle based on coordinate system
    if strcmp(ds.coord_sys, "earth")
        % Convert "deg CCW from East" to "deg CW from North" [0, 360]
        angle_compass = -(angle_deg - 90);
        relative_to = ds.vel.coords.dir{N_idx};  % North component
    else
        % Switch to clockwise and from [-180, 180] to [0, 360]
        angle_compass = -angle_deg;
        angle_compass(angle_compass < 0) = angle_compass(angle_compass < 0) + 360;
        relative_to = ds.vel.coords.dir{E_idx};  % East/streamwise component
    end

    % Ensure angles are in [0, 360] range
    angle_compass = mod(angle_compass, 360);
    % Set direction to NaN where velocity is zero
    angle_compass(U_mag == 0) = NaN;

    ds_out.U_dir = struct();
    % Create direction output
    ds_out.U_dir.data = angle_compass;
    ds_out.U_dir.dims = ds.vel.dims(1:end ~= dir_idx);
    ds_out.U_dir.coords = ds.vel.coords;
    ds_out.U_dir.units = "degrees clockwise from " + relative_to;
    ds_out.U_dir.standard_name = "sea_water_to_direction";
    ds_out.U_dir.long_name = "Water To Direction";
end
