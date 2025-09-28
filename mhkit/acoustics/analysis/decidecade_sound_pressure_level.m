function mspl = decidecade_sound_pressure_level(spsd, fmin, fmax)
% 
% Calculates the sound pressure level in decidecade bands directly
% from the mean square sound pressure spectral density (SPSD).
% 
% Parameters
% ----------
% spsd: struct
%     Mean square sound pressure spectral density in [Pa^2/Hz].
% fmin: int
%     Lower frequency band limit (lower limit of the hydrophone).
%     Default: 10 Hz
% fmax: int
%     Upper frequency band limit (Nyquist frequency).
%     Default: 100000 Hz
% 
% Returns
% -------
% mspl : struct
%     Sound pressure level [dB re 1 uPa] indexed by time and decidecade bands


arguments (Input)
    spsd struct
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
end

arguments (Output)
    mspl struct
end

fmax = fmax_warning(spsd.fs/2, fmax);

octave = 10;
base = 10;

mspl = band_sound_pressure_level(spsd, octave, base, fmin, fmax);
mspl.units = 'dB re 1 uPa';
mspl.name = 'Decidecade Sound Pressure Level';

end