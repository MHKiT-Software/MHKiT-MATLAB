function [slope, intercept] = check_turbulence_cascade_slope(psd, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the slope of the PSD within a frequency range to verify 
% isotropic turbulence cascade follows f^(-5/3) scaling
%
% This function calculates the slope of the power spectral density (PSD)
% within a given frequency range using log-log linear regression. The
% purpose is to check that the region containing the isotropic turbulence
% cascade decreases at the expected rate of f^(-5/3).
%
% Parameters
% ------------
%   psd : structure  
%       Power spectral density structure with:
%       psd.data : PSD values [time x freq] or [range x time x freq]
%       psd.coords.freq : frequency vector [Hz or rad/s]
%   freq_range : array, optional (name-value)
%       [min_freq, max_freq] range for cascade analysis (default: [0.2, 0.4])
%
% Returns
% ---------
%   slope : double array
%       Slope coefficients from log-log regression (should be ~-5/3)
%       Size matches PSD dimensions excluding frequency
%   intercept : double array  
%       Intercept coefficients from log-log regression
%       Size matches PSD dimensions excluding frequency
%
% Key Equations
% -------------
% 1. Isotropic turbulence cascade:
%    S(f) = α ε^(2/3) f^(-5/3) + N
%
% 2. Log-log linear regression:
%    log10(S) = m*log10(f) + b
%    where m should ≈ -5/3 for isotropic cascade
%
% Notes
% -----
% The slope is calculated using linear regression on log-transformed data:
% log10(PSD) vs log10(frequency). For ideal isotropic turbulence cascade,
% the slope should be -5/3 ≈ -1.667.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    psd struct
    options.freq_range (1,2) {mustBeNumeric} = [0.2, 0.4]
end

% Extract parameters
freq_range = options.freq_range;

% Validate input structure
if ~isfield(psd, 'data') || ~isfield(psd, 'coords') || ~isfield(psd.coords, 'freq')
    error('check_turbulence_cascade_slope:InvalidPSD', ...
        'psd must have .data and .coords.freq fields');
end

psd_data = psd.data;
freq = psd.coords.freq;

% Validate frequency range
if freq_range(1) >= freq_range(2)
    error('check_turbulence_cascade_slope:InvalidRange', ...
        'freq_range(1) must be less than freq_range(2)');
end

% Find frequency indices within range
idx = (freq > freq_range(1)) & (freq < freq_range(2));

if sum(idx) < 3
    error('check_turbulence_cascade_slope:InsufficientData', ...
        'Insufficient frequency points in range [%.3f, %.3f]. Need at least 3 points.', ...
        freq_range(1), freq_range(2));
end

% Extract frequencies and PSD values within range
freq_subset = freq(idx);
x = log10(freq_subset);  % log10(frequency)

% Handle different PSD dimensions
psd_dims = size(psd_data);
n_dims = length(psd_dims);

if n_dims == 2
    % 2D: [time x freq] 
    psd_subset = psd_data(:, idx);  % [time x freq_subset]
    y = log10(psd_subset);          % log10(PSD)
    
    % Calculate means for regression
    x_bar = mean(x);                % scalar
    y_bar = mean(y, 2);             % [time x 1]
    
    % Linear regression: y = mx + b
    % m = sum((x - x_bar) * (y - y_bar)) / sum((x - x_bar)^2)
    x_centered = x - x_bar;         % [freq_subset x 1]
    y_centered = y - y_bar;         % [time x freq_subset]
    
    % Transpose x_centered for proper broadcasting with 2D y_centered
    x_centered_t = x_centered';     % [1 x freq_subset]
    
    numerator = sum(x_centered_t .* y_centered, 2);    % [time x 1]
    denominator = sum(x_centered.^2);                  % scalar
    
    slope = numerator / denominator;                   % [time x 1]
    intercept = y_bar - slope * x_bar;                 % [time x 1]
    
elseif n_dims == 3
    % 3D: [range x time x freq]
    psd_subset = psd_data(:, :, idx);  % [range x time x freq_subset]
    y = log10(psd_subset);             % log10(PSD)
    
    % Calculate means for regression  
    x_bar = mean(x);                   % scalar
    y_bar = mean(y, 3);                % [range x time]
    
    % Linear regression along frequency dimension
    x_centered = x - x_bar;            % [freq_subset x 1]
    y_centered = y - y_bar;            % [range x time x freq_subset]
    
    % Reshape x_centered for proper broadcasting with 3D y_centered
    x_centered_3d = reshape(x_centered, [1, 1, length(x_centered)]);  % [1 x 1 x freq_subset]
    
    numerator = sum(x_centered_3d .* y_centered, 3);   % [range x time]
    denominator = sum(x_centered.^2);                  % scalar
    
    slope = numerator / denominator;                   % [range x time]
    intercept = y_bar - slope * x_bar;                 % [range x time]
    
else
    error('check_turbulence_cascade_slope:UnsupportedDimensions', ...
        'PSD data must be 2D [time x freq] or 3D [range x time x freq]');
end

% Ensure outputs are double precision
slope = double(slope);
intercept = double(intercept);

end