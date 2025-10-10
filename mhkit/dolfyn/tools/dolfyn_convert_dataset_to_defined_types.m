function ds = dolfyn_convert_dataset_to_defined_types(ds)
% Convert dataset fields to their defined types from instrument definitions
%
% This function converts a dolfyn dataset to use the correct data types
% as specified in the instrument definition files, ensuring consistency
% with NetCDF files and Python processing.
%
% Parameters
% ----------
% ds : struct
%     Dataset structure from dolfyn_read
%
% Returns
% -------
% ds : struct
%     Dataset with fields converted to defined types

    % Determine instrument type and get type definitions
    type_map = get_instrument_type_definitions(ds);

    % Apply type conversions
    ds = apply_type_conversions(ds, type_map);
end

function type_map = get_instrument_type_definitions(ds)
    % Dispatch to appropriate instrument definitions

    % Try to determine instrument type from dataset
    inst_type = detect_instrument_type(ds);

    switch lower(inst_type)
        case 'nortek'
            type_map = get_nortek_type_map();
        case 'signature'
            type_map = get_signature_type_map();
        case 'rdi'
            type_map = get_rdi_type_map();
        otherwise
            warning('Unknown instrument type: %s. Using signature defaults.', inst_type);
            type_map = get_signature_type_map();
    end
end

function inst_type = detect_instrument_type(ds)
    % Detect instrument type from dataset attributes

    % Check for signature-specific fields
    if isfield(ds, 'burst_config') || isfield(ds, 'attrs') && isfield(ds.attrs, 'burst_config')
        inst_type = 'signature';
        return;
    end

    % Check for RDI-specific fields
    if isfield(ds, 'attrs') && isfield(ds.attrs, 'sourceprog') && contains(lower(ds.attrs.sourceprog), 'rdi')
        inst_type = 'rdi';
        return;
    end

    % Check for Nortek Vector/AWAC fields
    if isfield(ds, 'attrs') && isfield(ds.attrs, 'inst_model') && contains(lower(ds.attrs.inst_model), 'vector')
        inst_type = 'nortek';
        return;
    end

    % Default to signature (most common case for our test files)
    inst_type = 'signature';
end

function type_map = get_signature_type_map()
    % Get type definitions for Signature instruments
    % Based on format codes from read_signature.m

    type_map = containers.Map();

    % Main data arrays (from comparison results)
    type_map('corr.data') = 'uint8';       % format 'B'
    type_map('xmit_energy.data') = 'uint16';   % format 'H'
    type_map('error.data') = 'uint16';     % format 'H'
    type_map('status.data') = 'uint32';    % format 'I'
    type_map('ensemble.data') = 'uint32';  % format 'I'
    type_map('error_b5.data') = 'uint16';  % format 'H'
    type_map('status_b5.data') = 'uint32'; % format 'I'
    type_map('ensemble_b5.data') = 'uint32'; % format 'I'
    type_map('corr_b5.data') = 'uint8';    % format 'B'
    type_map('xmit_energy_b5.data') = 'uint16'; % format 'H'
    type_map('nsamp_altraw.data') = 'uint32';   % format 'I'

    % Coordinate arrays (int32 from Python control files)
    type_map('coords.beam') = 'int32';
    type_map('coords.dir') = 'int32';
    type_map('coords.x1') = 'int32';
    type_map('coords.x2') = 'int32';
    type_map('vel.coords.dir') = 'int32';
    type_map('amp.coords.beam') = 'int32';
    type_map('corr.coords.beam') = 'int32';
    type_map('beam2inst_orientmat.coords.x1') = 'int32';
    type_map('beam2inst_orientmat.coords.x2') = 'int32';
end

function type_map = get_nortek_type_map()
    % Get type definitions from nortek_defs function - using existing source of truth

    type_map = containers.Map();

    % Get nortek definitions
    [vec_data, vec_sysdata, awac_profile] = nortek_defs();

    % Extract types from vec_data
    fields = fieldnames(vec_data);
    for i = 1:length(fields)
        field = fields{i};
        if isfield(vec_data.(field), 'dtype')
            type_map(sprintf('%s.data', field)) = char(vec_data.(field).dtype);
        end
    end

    % Extract types from vec_sysdata
    fields = fieldnames(vec_sysdata);
    for i = 1:length(fields)
        field = fields{i};
        if isfield(vec_sysdata.(field), 'dtype')
            type_map(sprintf('%s.data', field)) = char(vec_sysdata.(field).dtype);
        end
    end

    % Extract types from awac_profile
    fields = fieldnames(awac_profile);
    for i = 1:length(fields)
        field = fields{i};
        if isfield(awac_profile.(field), 'dtype')
            type_map(sprintf('%s.data', field)) = char(awac_profile.(field).dtype);
        end
    end
end

function type_map = get_rdi_type_map()
    % Get type definitions for RDI instruments from read_rdi.m definitions

    type_map = containers.Map();

    % Extract types from read_rdi.m (these are defined in that file)
    type_map('vel.data') = 'float32';      % per read_rdi.m line 173
    type_map('amp.data') = 'uint8';        % per read_rdi.m line 176
    type_map('corr.data') = 'uint8';       % per read_rdi.m line 180
    type_map('prcnt_gd.data') = 'uint8';   % per read_rdi.m line 184
    type_map('status.data') = 'float32';   % per read_rdi.m line 188
    type_map('number.data') = 'uint32';    % per read_rdi.m line 101
    type_map('rtc.data') = 'uint16';       % per read_rdi.m line 104
    type_map('adc.data') = 'uint8';        % per read_rdi.m line 152
end

function ds = apply_type_conversions(ds, type_map)
    % Apply type conversions to dataset based on type map

    field_paths = keys(type_map);

    for i = 1:length(field_paths)
        field_path = field_paths{i};
        target_type = type_map(field_path);

        % Navigate to the field and convert its type
        try
            current_value = get_nested_field(ds, field_path);
            if ~isempty(current_value)
                converted_value = convert_to_type(current_value, target_type);
                ds = set_nested_field(ds, field_path, converted_value);
            end
        catch ME
            % Field doesn't exist or can't be converted - skip it
            continue;
        end
    end
end

function value = get_nested_field(ds, field_path)
    % Get value from nested field path like 'corr.data'
    parts = strsplit(field_path, '.');
    value = ds;
    for i = 1:length(parts)
        if isfield(value, parts{i})
            value = value.(parts{i});
        else
            value = [];
            return;
        end
    end
end

function ds = set_nested_field(ds, field_path, value)
    % Set value in nested field path like 'corr.data'
    parts = strsplit(field_path, '.');

    if length(parts) == 1
        % Top-level field
        ds.(parts{1}) = value;
    elseif length(parts) == 2
        % Two-level field like 'corr.data'
        if isfield(ds, parts{1})
            ds.(parts{1}).(parts{2}) = value;
        end
    elseif length(parts) == 3
        % Three-level field like 'vel.coords.dir'
        if isfield(ds, parts{1}) && isfield(ds.(parts{1}), parts{2})
            ds.(parts{1}).(parts{2}).(parts{3}) = value;
        end
    else
        % For deeper nesting, use recursive approach
        if isfield(ds, parts{1})
            sub_path = strjoin(parts(2:end), '.');
            ds.(parts{1}) = set_nested_field(ds.(parts{1}), sub_path, value);
        end
    end
end

function converted_value = convert_to_type(value, target_type)
    % Convert value to target type

    switch lower(target_type)
        case 'uint8'
            converted_value = uint8(value);
        case 'uint16'
            converted_value = uint16(value);
        case 'uint32'
            converted_value = uint32(value);
        case 'int8'
            converted_value = int8(value);
        case 'int16'
            converted_value = int16(value);
        case 'int32'
            converted_value = int32(value);
        case 'single'
            converted_value = single(value);
        case 'float32'
            converted_value = single(value);
        case 'double'
            converted_value = double(value);
        case 'float64'
            converted_value = double(value);
        otherwise
            warning('Unknown target type: %s', target_type);
            converted_value = value;
    end
end


function ds = apply_standardization_fixes(ds)
    % Apply standardization fixes to match Python control file format
    % This function handles common differences between MATLAB and Python outputs

    % Fix units inconsistencies (Python uses "1" for dimensionless, MATLAB varies)
    ds = fix_units_standardization(ds);

    % Convert logical flags to uint8 for consistency with Python
    ds = convert_logical_flags_to_uint8(ds);

    % Fix dimension orientation (transpose issues)
    ds = fix_dimension_orientations(ds);

    % Fix samp_altraw data format
    ds = fix_samp_altraw_format(ds);

    % Fix burst_config struct vs char format
    ds = fix_burst_config_format(ds);

    % Add missing fields for Python control compatibility
    ds = add_missing_fields_for_python_control(ds);
end


function ds = fix_units_standardization(ds)
    % Standardize units to match Python control files
    % Python uses "1" for dimensionless quantities, MATLAB often uses specific units

    % List of fields that should have units "1" (dimensionless)
    dimensionless_fields = {'amp', 'amp_b5', 'beam2inst_orientmat', 'orientmat'};

    for i = 1:length(dimensionless_fields)
        field_name = dimensionless_fields{i};
        if isfield(ds, field_name) && isfield(ds.(field_name), 'units')
            % Change units to "1" for consistency with Python
            ds.(field_name).units = '1';
        end
    end
end


function ds = convert_logical_flags_to_uint8(ds)
    % Convert logical flag fields to uint8 for consistency with Python control files
    % Python stores boolean flags as uint8 (0/1), MATLAB often uses logical

    % List of flag fields that should be uint8 instead of logical
    flag_fields = {'low_volt_skip', 'active_config', 'telemetry_data', 'boost_running'};

    for i = 1:length(flag_fields)
        field_name = flag_fields{i};
        if isfield(ds, field_name) && isfield(ds.(field_name), 'data')
            if islogical(ds.(field_name).data)
                % Convert logical to uint8
                ds.(field_name).data = uint8(ds.(field_name).data);
            end
        end
    end
end


function ds = fix_dimension_orientations(ds)
    % Fix dimension orientation differences (transpose issues)
    % Python control files sometimes have different dimension orientations

    % Fields that commonly have dimension orientation issues
    transpose_fields = {'samp_altraw', 'beam2inst_orientmat'};

    for i = 1:length(transpose_fields)
        field_name = transpose_fields{i};
        if isfield(ds, field_name) && isfield(ds.(field_name), 'dims')
            % Check if dims is a cell array and ensure correct orientation
            if iscell(ds.(field_name).dims)
                % Convert to row vector format to match Python control files
                if size(ds.(field_name).dims, 1) > size(ds.(field_name).dims, 2)
                    ds.(field_name).dims = ds.(field_name).dims';
                end
            end
        end
    end
end


function ds = fix_samp_altraw_format(ds)
    % Fix samp_altraw data format and dimensions
    % Python control: double [1 1 1490], MATLAB: single [1490 1]

    if isfield(ds, 'samp_altraw') && isfield(ds.samp_altraw, 'data')
        % Convert from single to double for consistency
        if isa(ds.samp_altraw.data, 'single')
            ds.samp_altraw.data = double(ds.samp_altraw.data);
        end

        % Fix dimension shape - Python expects [1 1 N] format
        current_size = size(ds.samp_altraw.data);
        if length(current_size) == 2 && current_size(1) > 1 && current_size(2) == 1
            % Reshape from [N 1] to [1 1 N] to match Python control
            ds.samp_altraw.data = reshape(ds.samp_altraw.data, [1, 1, current_size(1)]);
        end
    end

    % Add missing coordinate and unit fields for altraw data if they don't exist
    altraw_fields = {'nsamp_altraw', 'dsamp_altraw', 'samp_altraw'};
    for i = 1:length(altraw_fields)
        field_name = altraw_fields{i};
        if isfield(ds, field_name)
            % Add coords.time if missing
            if ~isfield(ds.(field_name), 'coords')
                ds.(field_name).coords = struct();
            end
            if ~isfield(ds.(field_name).coords, 'time')
                % Use main time coordinate if available
                if isfield(ds, 'coords') && isfield(ds.coords, 'time')
                    ds.(field_name).coords.time = ds.coords.time;
                end
            end

            % Add units if missing
            if ~isfield(ds.(field_name), 'units')
                ds.(field_name).units = 'm';  % Default units for altimeter data
            end
        end
    end

    % Add n_altraw field if missing (derived from samp_altraw size)
    if isfield(ds, 'samp_altraw') && ~isfield(ds, 'n_altraw')
        if isfield(ds.samp_altraw, 'data')
            altraw_size = size(ds.samp_altraw.data);
            % n_altraw is typically the length dimension
            ds.n_altraw = max(altraw_size);
        end
    end
end


function ds = fix_burst_config_format(ds)
    % Fix burst_config format difference
    % Python expects burst_config as char, MATLAB creates it as struct

    if isfield(ds, 'attrs')
        % Find all burst_config fields in attrs
        attr_fields = fieldnames(ds.attrs);
        burst_config_fields = attr_fields(contains(attr_fields, 'burst_config'));

        for i = 1:length(burst_config_fields)
            field_name = burst_config_fields{i};
            if isstruct(ds.attrs.(field_name))
                % Convert struct to a char representation for Python compatibility
                % This is a simple approach - convert to JSON-like string
                config_struct = ds.attrs.(field_name);
                config_str = struct_to_config_string(config_struct);
                ds.attrs.(field_name) = config_str;
            end
        end
    end
end


function config_str = struct_to_config_string(config_struct)
    % Convert burst config struct to string representation
    % Python expects a string format for burst_config

    field_names = fieldnames(config_struct);
    config_parts = {};

    for i = 1:length(field_names)
        field_name = field_names{i};
        field_value = config_struct.(field_name);

        if islogical(field_value)
            % Convert logical to string
            if field_value
                value_str = 'true';
            else
                value_str = 'false';
            end
        elseif isnumeric(field_value)
            % Convert numeric to string
            value_str = num2str(field_value);
        else
            % Already a string
            value_str = char(field_value);
        end

        config_parts{end+1} = sprintf('%s=%s', field_name, value_str);
    end

    % Join all parts with semicolons (or another appropriate separator)
    config_str = strjoin(config_parts, ';');
end