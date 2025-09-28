function validate_spsd_struct(spsd_struct, function_name, options)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Validates SPSD (Sound Pressure Spectral Density) structure format and content
%
% Parameters
% ------------
%   spsd_struct: struct
%       SPSD structure to validate
%   function_name: string
%       Name of calling function for error messages
%   options: struct (optional)
%       require_positive : logical - Check data > 0 for log operations (default: false)
%       check_dimensions : logical - Check dimensional consistency (default: true)
%       required_fields : cell array - Additional required fields beyond 'data' (default: {})
%
% Returns
% ---------
%   No return value - throws error if validation fails
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd_struct struct
    function_name {mustBeText}
    options.require_positive {mustBeNumericOrLogical} = false
    options.check_dimensions {mustBeNumericOrLogical} = true
    options.required_fields cell = {}
end

% Always require 'data' field
base_required_fields = {'data'};
all_required_fields = [base_required_fields, options.required_fields];

% Validate required fields exist
for i = 1:length(all_required_fields)
    field_name = all_required_fields{i};
    if ~isfield(spsd_struct, field_name)
        error(sprintf('MHKiT:acoustics:%s:MissingField', function_name), ...
            'spsd structure must contain %s field', field_name);
    end
end

% Validate data field
if ~isnumeric(spsd_struct.data) || ~all(isfinite(spsd_struct.data), 'all')
    error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
        'spsd.data must be numeric and finite');
end

% Check positive values if required (for log operations)
if options.require_positive && any(spsd_struct.data <= 0, 'all')
    error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
        'spsd.data must contain positive values for log calculation');
end

% Validate freq field if present
if isfield(spsd_struct, 'freq')
    if ~isnumeric(spsd_struct.freq) || ~isvector(spsd_struct.freq) || any(spsd_struct.freq <= 0)
        error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
            'spsd.freq must be a positive numeric vector');
    end
end

% Validate time field if present
if isfield(spsd_struct, 'time')
    if ~(isnumeric(spsd_struct.time) || isdatetime(spsd_struct.time))
        error(sprintf('MHKiT:acoustics:%s:InvalidInput', function_name), ...
            'spsd.time must be numeric or datetime');
    end
end

% Check dimensional consistency if requested
if options.check_dimensions
    if isfield(spsd_struct, 'freq') && size(spsd_struct.data, 1) ~= length(spsd_struct.freq)
        error(sprintf('MHKiT:acoustics:%s:DimensionMismatch', function_name), ...
            'spsd.data rows (%d) must match spsd.freq length (%d)', ...
            size(spsd_struct.data, 1), length(spsd_struct.freq));
    end

    if isfield(spsd_struct, 'time') && size(spsd_struct.data, 2) ~= length(spsd_struct.time)
        error(sprintf('MHKiT:acoustics:%s:DimensionMismatch', function_name), ...
            'spsd.data columns (%d) must match spsd.time length (%d)', ...
            size(spsd_struct.data, 2), length(spsd_struct.time));
    end
end

end