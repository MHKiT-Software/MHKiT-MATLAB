function spl = sound_pressure_level(spsd, fmin, fmax)
% 
% Calculates the sound pressure level (SPL) in a specified frequency band
% from the mean square sound pressure spectral density (SPSD).
% 
% Parameters
% ----------
% spsd: struct
%     Mean square sound pressure spectral density in [Pa^2/Hz]
% fmin: integer
%     Lower frequency band limit (lower limit of the hydrophone). Default: 10 Hz
% fmax: integer
%     Upper frequency band limit (Nyquist frequency). Default: 100000 Hz
% 
% Returns
% -------
% spl: struct
%     Sound pressure level [dB re 1 uPa]

arguments (Input)
    spsd struct
    fmin {mustBeInteger} = 10
    fmax {mustBeInteger} = 100000
end

arguments (Output)
    spl struct
end

% type checks
if fmin <= 0
    error('fmin must be positive');
end
if fmax <= fmin
    error('fmax must be greater than fmin');
end
fmax = fmax_warning(spsd.fs/2, fmax);

% Reference value of sound pressure
reference = 1e-12; % Pa^2, = 1 uPa^2

% Select frequency band
idx = spsd.freq >= fmin & spsd.freq <= fmax;
band = spsd.data(idx, :);
freqs = spsd.freq(idx);

% Integrate mean square sound pressure over frequency band
pressure_squared = trapz(freqs, band);

% Mean square sound pressure level
mspl = 10 * log10(pressure_squared / reference);

% Output
spl.data = mspl;
spl.time = spsd.time;
spl.name = 'Sound Pressure Level';
spl.units = 'db re 1 uPa';

end
