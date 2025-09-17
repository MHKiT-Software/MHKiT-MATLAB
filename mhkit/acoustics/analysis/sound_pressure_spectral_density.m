function spsd = sound_pressure_spectral_density(data, fs, bin_length)
% 
% Calculates the sound pressure spectral density (SPSD) from audio
% samples split into FFTs with a specified bin length in seconds,
% using Hanning windowing with 50% overlap.
% 
% Parameters
% ----------
% data : struct
%     Sound pressure in [Pa] or Voltage [V]
% fs : double
%     Data collection sampling rate [Hz]
% bin_length: double
%     Length of time in seconds to create FFTs. Default: 1.
% 
% Returns
% -------
% spsd: struct
%     Spectral density [Pa^2/Hz] indexed by time and frequency

arguments (Input)
    data struct
    fs {mustBeNumeric}
    bin_length {mustBeNumeric} = 1
end

arguments (Output)
    spsd struct
end

pressure = data.data;

% cut into multiple 1-second samples, calculate PSD of each sample
win  = bin_length*fs;  % window length of each time series
step = 0.5*win; % overlap betewen each window
ns   = floor((length(pressure)-win)/step)+1; % number of time series samples
nfft = win; % number of fft points
df   = fs/nfft;  % frequency resolution
nfreq= ceil((nfft)/2);  % Next highest power of 2 greater than length(x).
freq = (1:nfreq)*2*fs/2/nfft; % frequency of spectrum
Pf2  = zeros(nfreq,ns);  % mean-squared sound pressure spectral density

for i=1:ns
    sample   = pressure((i-1)*step+1:(i-1)*step+win);
    sample   = sample-mean(sample); 
    t_power  = sum(sample.^2)/length(sample); % mean squred sound pressure; power in time domain
    [f,spec] = fft_hann(fs,sample, nfft);  % spectrum
    psd      = spec.*conj(spec)/df/2;  % PSD
    f_power  = sum(psd)*df; % power in frequency domain
    psd_adj  = psd*t_power/f_power; % adjust the amplitude of PSD according to Parseval's theorem
    Pf2(:,i) = psd_adj; %mean-squared sound pressure spectral density
end

spsd.data = Pf2;
spsd.time = linspace(data.time(1), data.time(end), ns);
spsd.freq = freq';
spsd.name = "Mean Square Sound Pressure Spectral Density";
spsd.units = strcat(data.units, "^2/Hz");
spsd.fs = fs;
spsd.nbin = strcat(string(bin_length)," s");
spsd.overlap = "50%";
spsd.nfft = nfft;

end