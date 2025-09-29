function [f, spec] = fft_hann(fs, x, nfft)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply Fast-Fourier-Transform to a time series using a Hanning window
%
% Parameters
% ------------
%   fs: double
%       Sample rate of data [Hz]
%   x: double array
%       Array of time series data
%   nfft: double
%       Number of elements in the FFT
%
% Returns
% ---------
%   f: double array
%       Frequency array vector [Hz]
%   spec: complex double array
%       Frequency spectra resulting from the FFT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    fs {mustBeNumeric, mustBePositive}
    x {mustBeNumeric}
    nfft {mustBeNumeric, mustBePositive, mustBeInteger}
end

arguments (Output)
    f {mustBeNumeric}
    spec {mustBeNumeric}
end

% Pre-process the input signal
x0 = x - mean(x);                                    % Demean the time series
N = length(x);                                       % Get length of time series

% Generate and apply Hanning window
w = 0.5 * (1 - cos(2*pi*(0:1:N-1) / (N-1)));         % Generate Hanning window coefficients
xx = x0 .* w';                                       % Apply Hanning window

% Calculate FFT parameters
Fn = fs / 2;                                         % Nyquist frequency
NumFFT = nfft;                                       % Number of FFT points
NumUniquePts = ceil(NumFFT / 2);                     % Number of unique FFT points

% Compute FFT
TempFFT = fft(xx, NumFFT);                           % Take FFT, padding with zeros if needed
spec = TempFFT(1:NumUniquePts);                      % Keep only positive frequencies

% Adjust spectrum for single-sided representation
spec = spec * 2;                                     % Account for negative frequencies
spec(1) = spec(1) / 2;                               % DC component should not be doubled
if mod(NumFFT, 2) == 0                               % For even NumFFT
    spec(end) = spec(end) / 2;                       % Nyquist frequency should not be doubled
end

% Normalize spectrum and generate frequency vector
spec = spec / length(xx);                           % Normalize by signal length
f = (1:NumUniquePts) * 2 * Fn / NumFFT;             % Generate frequency vector

end
