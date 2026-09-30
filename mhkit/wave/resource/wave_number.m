function k = wave_number(f, h, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave number from frequency and water depth
%
% Solves the linear dispersion relation (Eq 11 in IEC 62600-101 Ed. 2.0 en 2024)
%
% To compute wave number from angular frequency (w), convert w to f before
% using this function (f = w / (2*pi))
%
% Parameters
% ------------
% f : vector or scalar [Hz]
%   Frequency
% h : double [m]
%   Water depth
% rho : double [kg/m^3] (optional)
%   Water density. Default = 1025 kg/m^3
% g : double [m/s^2] (optional)
%   Gravitational acceleration. Default = 9.80665 m/s^2
%
% Returns
% ---------
% k : struct
%   k.values : vector [1/m]
%     Wave number
%   k.frequency : vector [Hz]
%     Frequency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    f {mustBeNumeric}
    h (1,1) {mustBeNumeric, mustBePositive}
    options.rho (1,1) {mustBeNumeric, mustBePositive} = 1025
    options.g (1,1) {mustBeNumeric, mustBePositive} = 9.80665
end

arguments (Output)
    k struct
end

g = options.g;

% Ensure f is a column vector for consistent processing
f_input = f;
f = f(:);

% Angular frequency
w = 2 * pi * f;
% note: = h*wa/sqrt(h*g/h)
xi = w / sqrt(g / h);
yi = xi.^2 ./ (1 - exp(-xi.^2.4908)).^0.4015;
k0 = yi / h;

% Solve dispersion relation: w^2 = g*k*tanh(k*h)
% Rearranged as: w^2 - g*k*tanh(k*h) = 0
k_values = k0;

% Set solver options to suppress output
fzero_options = optimset('Display', 'off', 'TolX', 1e-12);

% Only solve for points where initial guess isn't accurate enough
for i = 1:length(f)
    % Eq 11 in IEC 62600-101 Ed. 2.0 en 2024 using initial guess from Guo (2002)
    residual = w(i)^2 - g * k0(i) * tanh(k0(i) * h);
    if abs(residual) > 1e-9
        func = @(kk) w(i)^2 - g * kk * tanh(kk * h);
        try
            k_values(i) = fzero(func, k0(i), fzero_options);
        catch ME
            error('MHKiT:wave_number:SolverFailed', ...
                'Wave number solver failed for f=%g Hz: %s', f(i), ME.message);
        end
    end
end

% Reshape output to match input shape
if isrow(f_input)
    k_values = k_values';
end

% Build output structure
k = struct();
k.values = k_values;
k.frequency = f_input;

end
