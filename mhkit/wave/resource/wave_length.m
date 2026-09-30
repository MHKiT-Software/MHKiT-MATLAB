function l = wave_length(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave length from wave number
%
% Computes: l = 2*pi/k
%
% Parameters
% ------------
% k : numeric (scalar, vector, matrix, or struct)
%   Wave number [1/m]. If struct, uses k.values field.
%
% Returns
% ---------
% l : numeric [m]
%   Wave length. Same type/shape as input k (or k.values if struct).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    k
end

arguments (Output)
    l {mustBeNumeric}
end

% Handle struct input (k.values)
if isstruct(k)
    if ~isfield(k, 'values')
        error('MHKiT:wave_length:InvalidInput', ...
            'Structure input must have a values field');
    end
    k_values = k.values;
else
    k_values = k;
end

% Validate input
if ~isnumeric(k_values)
    error('MHKiT:wave_length:InvalidInput', ...
        'k must be numeric or a struct with numeric values field');
end

l = 2 * pi ./ k_values;

end
