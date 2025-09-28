function validate_spsdl_struct(spsdl_struct, function_name, options)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Validates SPSDL (Sound Pressure Spectral Density Level) structure format and content
%
% Parameters
% ------------
%   spsdl_struct: struct
%       SPSDL structure to validate
%   function_name: string
%       Name of calling function for error messages
%   options: struct (optional)
%       require_positive : logical - Check data > 0 for log operations (default: false)
%       check_dimensions : logical - Check dimensional consistency (default: true)
%       required_fields : cell array - Override default required fields (default: {'data', 'freq', 'time', 'name', 'units'})
%
% Returns
% ---------
%   No return value - throws error if validation fails
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsdl_struct struct
    function_name {mustBeText}
    options.require_positive {mustBeNumericOrLogical} = false
    options.check_dimensions {mustBeNumericOrLogical} = true
    options.required_fields cell = {'data', 'freq', 'time', 'name', 'units'}
end

% Validate required fields exist
for i = 1:length(options.required_fields)
    field_name = options.required_fields{i};
    if ~isfield(spsdl_struct, field_name)
        error(sprintf('MHKiT:acoustics:%s:MissingField', function_name), ...
            'spsdl structure must contain %s field', field_name);
    end
end

% Validate data field
if ~isnumeric(spsdl_struct.data) || ~all(isfinite(spsdl_struct.data), 'all')
    error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
        'spsdl.data must be numeric and finite');
end

% Check positive values if required (for log operations)
if options.require_positive && any(spsdl_struct.data <= 0, 'all')
    error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
        'spsdl.data must contain positive values for log calculation');
end

% Validate freq field if present
if isfield(spsdl_struct, 'freq')
    if ~isnumeric(spsdl_struct.freq) || ~isvector(spsdl_struct.freq) || any(spsdl_struct.freq <= 0)
        error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
            'spsdl.freq must be a positive numeric vector');
    end
end

% Validate time field if present
if isfield(spsdl_struct, 'time')
    if ~(isnumeric(spsdl_struct.time) || isdatetime(spsdl_struct.time))
        error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
            'spsdl.time must be numeric or datetime');
    end
end

% Check dimensional consistency if requested
if options.check_dimensions
    if isfield(spsdl_struct, 'freq') && size(spsdl_struct.data, 1) ~= length(spsdl_struct.freq)
        error(sprintf('MHKiT:acoustics:%s:DimensionMismatch', function_name), ...
            'spsdl.data rows (%d) must match spsdl.freq length (%d)', ...
            size(spsdl_struct.data, 1), length(spsdl_struct.freq));
    end

    if isfield(spsdl_struct, 'time') && size(spsdl_struct.data, 2) ~= length(spsdl_struct.time)
        error(sprintf('MHKiT:acoustics:%s:DimensionMismatch', function_name), ...
            'spsdl.data columns (%d) must match spsdl.time length (%d)', ...
            size(spsdl_struct.data, 2), length(spsdl_struct.time));
    end
end

end