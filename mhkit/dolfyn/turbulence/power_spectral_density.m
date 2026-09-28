function ds = power_spectral_density(ds, vel_data, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate power spectral density from ADCP velocity data using Welch's method.
%
% Parameters
% ------------
%   ds : structure
%       Dataset structure to add PSD results to
%   vel_data : double
%       Velocity time series data (1D array) [m/s]
%   freq_units : string, optional (name-value)
%       Output frequency units: 'Hz' (default) or 'rad/s'
%   n_fft : double, optional (name-value)
%       FFT window size (default: length(vel_data)/3) [samples]
%   field_name : string, optional (name-value)
%       Name of output field (default: 'psd')
%   fs : double, optional (name-value)
%       Sampling frequency (default: from ds.attrs.fs) [Hz]
%   window : string, optional (name-value)
%       Window function: 'hann' (default), 'hamming', 'rectwin'
%   overlap : double, optional (name-value)
%       Overlap fraction 0-1 (default: 0.5)
%
% Returns
% ---------
%   ds : structure
%       Input dataset with added PSD field containing:
%       .data : PSD values [m^2/s^2/Hz] or [m^2/s^2/(rad/s)]
%       .coords.frequency : frequency coordinates
%       .attrs.units : 'm2 s-2 Hz-1' or 'm2 s rad-1'
%       .attrs.n_fft : FFT window size
%       .attrs.long_name : 'Power Spectral Density'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    ds (:,:) struct
    vel_data (:,1) {mustBeNumeric}
    options.freq_units {mustBeTextScalar} = 'Hz'
    options.n_fft {mustBeNumeric, mustBePositive} = []
    options.field_name {mustBeTextScalar} = 'psd'
    options.fs {mustBeNumeric, mustBePositive} = []
    options.window {mustBeTextScalar} = 'hann'
    options.overlap (1,1) {mustBeNumeric, mustBeNonnegative, mustBeLessThan(options.overlap,1)} = 0.5
end

% Extract parameters
freq_units = options.freq_units;
n_fft = options.n_fft;
field_name = options.field_name;
fs = options.fs;
window_type = options.window;
overlap = options.overlap;

% Get sampling frequency and ensure it's double
if isempty(fs)
    if isfield(ds, 'attrs') && isfield(ds.attrs, 'fs')
        fs = double(ds.attrs.fs);
    else
        error('MHKiT:dolfyn:power_spectral_density:NoSamplingFreq', 'Sampling frequency not found. Provide fs parameter or ds.attrs.fs');
    end
else
    fs = double(fs);
end

% Handle input - raw time series that needs to be binned like ds
vel_data = double(vel_data(:));
n_samples = length(vel_data);

% Get binning information from ds.attrs.n_bin (set by average_by_dimension)
if isfield(ds, 'attrs') && isfield(ds.attrs, 'n_bin')
    bin_size = double(ds.attrs.n_bin);
    n_time_bins = floor(n_samples / bin_size);
else
    error('MHKiT:dolfyn:power_spectral_density:NoNBin', 'ds.attrs.n_bin not found. Run average_by_dimension first.');
end

% Set default n_fft to match Python example: n_fft = ds.n_bin // 2
if isempty(n_fft)
    if isfield(ds, 'attrs') && isfield(ds.attrs, 'n_bin')
        n_fft = double(floor(ds.attrs.n_bin / 2));
    else
        n_fft = double(floor(n_samples / 3));  % Default fallback
    end
else
    n_fft = double(n_fft);
end

% Validate n_fft against bin size (not total data length)
if n_fft > bin_size
    warning('MHKiT:dolfyn:power_spectral_density:LargeNFFT', 'n_fft (%d) > bin_size (%d), using bin_size', n_fft, bin_size);
    n_fft = bin_size;
end

% Reshape data into time bins [bin_size x n_time_bins]
usable_samples = n_time_bins * bin_size;
vel_binned = reshape(vel_data(1:usable_samples), bin_size, n_time_bins);

% Calculate step size parameters
[step, n_segments, n_fft_used] = dolfyn_stepsize(bin_size, n_fft, [], []);

% Create frequency vector (will be computed by dolfyn_psd_1D, but we need it for output)
n_freq_end = floor(n_fft / 2.0 + 1); % Python: int(nfft / 2.0 + 1)
freq_indices = (2:n_freq_end)'; % Skip DC component (MATLAB index 1), start from index 2
freq = (freq_indices - 1) * fs / n_fft; % Convert to actual frequencies
n_freq = length(freq);

% Initialize output array
psd_vals = zeros(n_time_bins, n_freq);

% Process each time bin independently using dolfyn_psd_1D
for t = 1:n_time_bins
    vel_bin = vel_binned(:, t);  % Extract time bin
    
    try
        psd_bin = dolfyn_psd_1D(vel_bin, n_fft, fs, 'window', window_type, 'step', step);
        psd_vals(t, :) = psd_bin';  % Ensure row vector
    catch ME
        % Handle any errors (e.g., insufficient data)
        % warning('MHKiT:dolfyn:power_spectral_density:PSDBinError', 'Error computing PSD for time bin %d: %s', t, ME.message);
        psd_vals(t, :) = NaN;
    end
end

% Convert frequency units if requested
if strcmpi(freq_units, 'rad/s')
    freq = freq * 2 * pi;
    psd_vals = psd_vals / (2 * pi);  % Convert PSD units too
    units_str = 'm2 s-1 rad-1';
else
    units_str = 'm2 s-2 Hz-1';
end

ds.(field_name) = struct();
ds.(field_name).data = single(psd_vals);
ds.(field_name).dims = {'time', 'freq'};
ds.(field_name).coords = struct();
ds.(field_name).coords.freq = single(freq);
ds.(field_name).attrs = struct();
ds.(field_name).attrs.units = units_str;
ds.(field_name).attrs.n_fft = int32(n_fft);
ds.(field_name).attrs.long_name = 'Power Spectral Density';

end
