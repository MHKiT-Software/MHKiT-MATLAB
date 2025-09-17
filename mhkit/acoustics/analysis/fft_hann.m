function [f, spec] = fft_hann(fs, x, nfft)
% 
% Apply Fast-Fourier-Transform to a time series using a Hanning window.
% 
% Parameters
% -----------
% fs: int
%     Sample rate of data.
% x: array
%     Array of timeseries data.
% nfft: int
%     Number of elements in the FFT.
% 
% Returns
% -------
% f: array
%     Frequency array vector.
% spec: array
%     Frequency spectra resulting from the FFT.

arguments (Input)
    fs {mustBeNumeric}
    x {mustBeNumeric}
    nfft {mustBeNumeric}
end

arguments (Output)
    f {mustBeNumeric}
    spec {mustBeNumeric}
end

x0 = x-mean(x); % demean the time series
N = length(x); % get length of time series
w = 0.5 * (1-cos(2*pi*(0:1:N-1) / (N-1))); % generate hanning window
xx = x0.*w'; % apply hamming window
Fn = fs/2; % Nyquist frequency
NumFFT  = nfft;
TempFFT = fft(xx,NumFFT);                     % Take FFT, padding with zeros.
NumUniquePts = ceil((NumFFT)/2);
spec         = TempFFT(1:NumUniquePts);         % FFT is symmetric, throw away second half                 
spec         = spec*2;                          % Multiply by 2 to take into account the fact that we
                                                % threw out second half of TempFFT above
spec(1)           = spec(1)/2;            % Account for endpoint uniqueness
spec(length(spec))= spec(length(spec))/2;  % We know NumFFT is even
spec              = spec/length(xx);       % Scale the FFT so that it is not a function of the length of x.
f                 = (1:NumUniquePts)*2*Fn/NumFFT;

end