function spsd = sound_pressure_spectral_density(data, fs, bin_length, fft_length, rms)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound pressure spectral density (SPSD) from audio samples split into FFTs
% with a specified bin length in seconds, using Hanning windowing with 50% overlap.
%
% Parameters
% ------------
%   data: structure
%       data.data : Sound pressure [Pa] or Voltage [V]
%       data.time : Time vector [s] or datetime
%       data.units : Data units string
%   fs: double
%       Data collection sampling rate [Hz]
%   bin_length: double
%       Length of time in seconds to create FFTs. Default: 1.
%   fft_length: double
%       Length of FFT to use. If none, uses bin_length*fs. Default: None.        
%   rms: logical
%       If true, scales PSD by mean square of original signal. Default: true.
%
% Returns
% ---------
%   spsd: structure
%       spsd.data : Mean square sound pressure spectral density [Pa^2/Hz or V^2/Hz]
%       spsd.time : Time vector for spectral density bins [s] or datetime
%       spsd.freq : Frequency vector [Hz]
%       spsd.name : Description string
%       spsd.units : Units string [Pa^2/Hz or V^2/Hz]
%       spsd.fs : Sampling frequency [Hz]
%       spsd.bin_length : Bin length string [s]
%       spsd.overlap : Overlap percentage string
%       spsd.n_fft : Number of FFT points
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    data struct
    fs {mustBeNumeric}
    bin_length {mustBeNumeric} = 1
    fft_length {mustBeNumeric} = 0
    rms logical = true
end

arguments (Output)
    spsd struct
end

% Validate data structure has required fields
if ~isfield(data, 'data')
    error('MHKiT:acoustics:sound_pressure_spectral_density:MissingField', ...
        'data structure must contain data field');
end

if ~isfield(data, 'time')
    error('MHKiT:acoustics:sound_pressure_spectral_density:MissingField', ...
        'data structure must contain time field');
end

% Validate field types
if ~isnumeric(data.data)
    error('MHKiT:acoustics:sound_pressure_spectral_density:InvalidInput', ...
        'data.data must be numeric');
end

% Accept both numeric and datetime for time fields
if ~(isnumeric(data.time) || isdatetime(data.time))
    error('MHKiT:acoustics:sound_pressure_spectral_density:InvalidInput', ...
        'data.time must be numeric or datetime');
end

% Check dimensional consistency: data length must match time length
if length(data.data) ~= length(data.time)
    error('MHKiT:acoustics:sound_pressure_spectral_density:DimensionMismatch', ...
        'data.data length (%d) must match data.time length (%d)', ...
        length(data.data), length(data.time));
end

pressure = data.data;

% cut into multiple samples, calculate PSD of each sample
if fft_length == 0
    win  = bin_length*fs;  % window length of each time series
else
    win = fft_length;
end
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
    [f,spec] = fft_hann(fs,sample, nfft);  % spectrum
    psd      = spec.*conj(spec)/df/2;  % PSD
    if rms
        t_power  = sum(sample.^2)/length(sample); % mean squred sound pressure; power in time domain
        f_power  = sum(psd)*df; % power in frequency domain
        psd_adj  = psd*t_power/f_power; % adjust the amplitude of PSD according to Parseval's theorem
    else
        psd_adj  = psd;
    end
    Pf2(:,i) = psd_adj; %mean-squared sound pressure spectral density
end

spsd.data = Pf2;
spsd.time = linspace(data.time(1), data.time(end), ns);
spsd.freq = freq';
if rms
    spsd.name = "Mean Square Sound Pressure Spectral Density";
else
    spsd.name = "Sound Pressure Spectral Density";
end
spsd.units = strcat(data.units, "^2/Hz");
spsd.fs = fs;
spsd.bin_length = bin_length;
spsd.overlap = "50%";
spsd.n_fft = nfft;

end
