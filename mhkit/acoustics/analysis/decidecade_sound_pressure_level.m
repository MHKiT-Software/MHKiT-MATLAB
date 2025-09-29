function mspl = decidecade_sound_pressure_level(spsd, fmin, fmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound pressure level in decidecade bands directly
% from mean square sound pressure spectral density (SPSD).
%
% Parameters
% ------------
%   spsd: structure
%       spsd.data : Mean square sound pressure spectral density [Pa^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector [s]
%       spsd.fs : Sampling frequency [Hz]
%   fmin: double
%       Lower frequency band limit [Hz]. Default: 10
%   fmax: double
%       Upper frequency band limit [Hz]. Default: 100000
%
% Returns
% ---------
%   mspl: structure
%       mspl.data : Sound pressure level [dB re 1 uPa]
%       mspl.freq : Decidecade band center frequencies [Hz]
%       mspl.time : Time vector [s]
%       mspl.units : Units string
%       mspl.name : Descriptive name string
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


arguments (Input)
    spsd struct
    fmin {mustBeNumeric, mustBePositive} = 10
    fmax {mustBeNumeric, mustBePositive} = 100000
end

arguments (Output)
    mspl struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'decidecade_sound_pressure_level', ...
    'required_fields', {{'freq', 'time', 'fs'}});

% Validate frequency range
if fmin >= fmax
    error('MHKiT:acoustics:decidecade_sound_pressure_level:InvalidInput', ...
        'fmin must be less than fmax');
end

fmax = fmax_warning(spsd.fs/2, fmax);

octave = 10;
base = 10;

mspl = band_sound_pressure_level(spsd, octave, base, fmin, fmax);
mspl.units = 'dB re 1 uPa';
mspl.name = 'Decidecade Sound Pressure Level';

end
