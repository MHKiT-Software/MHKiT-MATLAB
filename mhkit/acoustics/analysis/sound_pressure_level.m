function spl = sound_pressure_level(spsd, fmin, fmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound pressure level (SPL) in a specified frequency band from
% the mean square sound pressure spectral density (SPSD)
%
% Parameters
% ------------
%   spsd: struct
%       spsd.data : Mean square sound pressure spectral density [Pa^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector
%       spsd.fs : Sampling frequency [Hz]
%   fmin: integer
%       Lower frequency band limit [Hz] (default: 10)
%   fmax: integer
%       Upper frequency band limit [Hz] (default: 100000)
%
% Returns
% ---------
%   spl: struct
%       spl.data : Sound pressure level [dB re 1 uPa]
%       spl.time : Time vector
%       spl.name : Data name identifier
%       spl.units : Data units string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    fmin {mustBeInteger} = 10
    fmax {mustBeInteger} = 100000
end

arguments (Output)
    spl struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'sound_pressure_level', ...
    'required_fields', {{'freq', 'time', 'fs'}});

% Additional parameter validation
if fmin <= 0
    error('MHKiT:acoustics:sound_pressure_level:InvalidInput', ...
        'fmin must be positive');
end
if fmax <= fmin
    error('MHKiT:acoustics:sound_pressure_level:InvalidInput', ...
        'fmax must be greater than fmin');
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
spl.units = 'dB re 1 uPa';

end
