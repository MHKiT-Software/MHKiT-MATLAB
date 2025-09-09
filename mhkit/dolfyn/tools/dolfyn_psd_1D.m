function psd = dolfyn_psd_1D(data, n_fft, fs, options)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Compute 1D power spectral density using Welch's method
%
%     Uses same algorithm as MHKiT-MHKiT-Python psd_1D function
%
% Parameters
% ------------
%   data : double
%       1D ADCP velocity time series data (column vector) [m/s]
%   n_fft : double
%       Number of points in the FFT window [samples]
%   fs : double
%       Sample rate [Hz]
%   window : char or double, optional (name-value)
%       Window function {'hann', 'hamming', 'rectwin'} or custom window vector (default: 'hann')
%   step : double, optional (name-value) 
%       Step size for segment overlap (default: auto-calculated) [samples]
%
% Returns
% ---------
%   psd : double
%       Power spectral density [m^2/s^2/Hz]
%       Size: [n_freq x 1] where n_freq = floor(n_fft/2)
%       Frequencies correspond to positive frequencies only (DC component excluded)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data (:,1) {mustBeNumeric}
    n_fft (1,1) {mustBeNumeric, mustBePositive}
    fs (1,1) {mustBeNumeric, mustBePositive}
    options.window = 'hann'
    options.step {mustBeNumeric} = []
end

% Extract optional parameters
window_type = options.window;
step = options.step;

% Ensure inputs are double and column vector
data = double(data(:));
n_fft = double(n_fft);
fs = double(fs);

% Get data length
l = length(data);

% Handle edge cases
if l < 2
    error('MHKiT:dolfyn_psd_1D:InsufficientData', 'Need at least 2 data points for PSD calculation');
end

% Check for mostly NaN data
valid_count = sum(~isnan(data));
if valid_count < 2
    error('MHKiT:dolfyn_psd_1D:InsufficientValidData', 'Need at least 2 valid (non-NaN) data points');
end

% Calculate loop step size and number of segments
[step, n_segments, n_fft] = dolfyn_stepsize(l, n_fft, [], step);

% Create window using extracted function
win = mhkit_get_window(window_type, n_fft);

% MHKiT-Python frequency indexing: slice(1, int(nfft / 2.0 + 1))
% This skips DC component (index 0) and goes to floor(nfft/2.0+1)-1
n_freq_end = floor(n_fft / 2.0 + 1);  % MHKiT-Python: int(nfft / 2.0 + 1)
freq_indices = (2:n_freq_end)';        % Skip DC (MATLAB index 1), start from 2
n_freq = length(freq_indices);

% Window normalization (MHKiT-Python: wght = 2.0 / (window**2).sum())
window_weight = 2.0 / sum(win.^2);

% Initialize PSD accumulator
psd_accumulator = zeros(n_freq, 1);

% MHKiT-Python algorithm: First segment outside loop, then loop for remaining
% Segment 0: MHKiT-Python s1 = fft(detrend_array(a[0:nfft]) * window)[fft_inds]
segment_start = 1;
segment_end = segment_start + n_fft - 1;
segment = data(segment_start:segment_end);

% Detrend segment
segment = mhkit_detrend_array(segment);

% Apply window and compute FFT
segment_windowed = segment .* win;
fft_segment = fft(segment_windowed, n_fft);

% Apply frequency indexing FIRST (like MHKiT-Python)
fft_selected = fft_segment(freq_indices);

% Compute power: MHKiT-Python pwr = np.abs(s1) ** 2
psd_accumulator = abs(fft_selected).^2;

% Loop for remaining segments (MHKiT-Python: if nens - 1)
if n_segments > 1
    % MHKiT-Python: for i in range(step, l - nfft + 1, step):
    for i = (step+1):step:(l - n_fft + 1)  % Convert to MATLAB 1-based indexing
        segment = data(i:(i + n_fft - 1));
        
        % Detrend segment
        segment = mhkit_detrend_array(segment);
        
        % Apply window and compute FFT
        segment_windowed = segment .* win;
        fft_segment = fft(segment_windowed, n_fft);
        
        % Apply frequency indexing FIRST
        fft_selected = fft_segment(freq_indices);
        
        % Accumulate power: MHKiT-Python pwr += np.abs(s1) ** 2
        psd_accumulator = psd_accumulator + abs(fft_selected).^2;
    end
end

% Apply final normalization (MHKiT-Python: pwr *= wght / nens / fs)
psd = psd_accumulator * window_weight / n_segments / fs;

end
