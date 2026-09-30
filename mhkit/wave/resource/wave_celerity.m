function Cg = wave_celerity(k, h, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave celerity (group velocity)
%
% Uses the formula from Eq 10 in IEC 62600-101 Ed. 2.0 en 2024:
% Cg = (pi * f / k) * (1 + (2*h*k) / sinh(2*h*k))
%
% For deep water (when depth_check=true and h/l > ratio), uses the
% simplified formula: Cg = pi * f / k
%
% Parameters
% ------------
% k : struct
%   Wave number structure:
%     k.values : vector [1/m]
%       Wave number values
%     k.frequency : vector [Hz]
%       Frequency
% h : double [m]
%   Water depth
% g : double [m/s^2] (optional)
%   Gravitational acceleration. Default = 9.80665 m/s^2
% depth_check : logical (optional)
%   If true, check depth regime and use deep water approximation
%   where applicable. Default = false
% ratio : double (optional)
%   Only applied if depth_check=true. If h/l > ratio,
%   water depth is set to deep. Default = 2
%
% Returns
% ---------
% Cg : struct
%   Cg.values : vector [m/s]
%     Wave celerity (group velocity)
%   Cg.frequency : vector [Hz]
%     Frequency
%   Cg.h : double [m]
%     Water depth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    k struct
    h (1,1) {mustBeNumeric, mustBePositive}
    options.g (1,1) {mustBeNumeric, mustBePositive} = 9.80665
    options.depth_check (1,1) logical = false
    options.ratio (1,1) {mustBeNumeric, mustBePositive} = 2
end

arguments (Output)
    Cg struct
end

% Validate input structure
if ~isfield(k, 'values') || ~isfield(k, 'frequency')
    error('MHKiT:wave_celerity:InvalidInput', ...
        'k must be a structure with values and frequency fields');
end

k_values = k.values(:);
f = k.frequency(:);

if options.depth_check
    % Calculate wavelength
    l = wave_length(k_values);

    % Get depth regime (true = deep water)
    dr = depth_regime(l, h, 'ratio', options.ratio);

    % Initialize output
    Cg_values = zeros(size(k_values));

    % Deep water approximation for deep frequencies
    if any(dr)
        Cg_values(dr) = pi * f(dr) ./ k_values(dr);
    end

    % Full formula for shallow/intermediate frequencies
    if any(~dr)
        sf = f(~dr);
        sk = k_values(~dr);
        Cg_values(~dr) = (pi * sf ./ sk) .* (1 + (2 * h * sk) ./ sinh(2 * h * sk));
    end
else
    % Eq 10 in IEC 62600-101 Ed. 2.0 en 2024
    Cg_values = (pi * f ./ k_values) .* (1 + (2 * h * k_values) ./ sinh(2 * h * k_values));
end

% Handle shape to match input
if isrow(k.values)
    Cg_values = Cg_values';
end

% Build output structure
Cg = struct();
Cg.values = Cg_values;
Cg.frequency = k.frequency;
Cg.h = h;

end
