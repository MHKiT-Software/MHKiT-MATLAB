function mspl = third_octave_sound_pressure_level(spsd, fmin, fmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound pressure level in third octave bands directly
% from the mean square sound pressure spectral density (SPSD).
%
% Parameters
% ------------
%     spsd : structure
%         Mean square sound pressure spectral density structure containing:
%         spsd.S : Mean square sound pressure spectral density [Pa^2/Hz]
%         spsd.f : Frequency vector [Hz]
%         spsd.fs : Sampling frequency [Hz]
%     fmin : double
%         Lower frequency band limit (typically hydrophone lower limit) [Hz]
%     fmax : double
%         Upper frequency band limit (typically Nyquist frequency) [Hz]
%
% Returns
% ---------
%     mspl : structure
%         Third octave band sound pressure level structure containing:
%         mspl.spl : Sound pressure level [dB re 1 μPa]
%         mspl.f : Center frequencies of third octave bands [Hz]
%         mspl.units : Units string 'dB re 1 uPa'
%         mspl.name : Descriptive name 'Third Octave Sound Pressure Level'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    fmin {mustBeNumeric, mustBePositive, mustBeFinite} = 10
    fmax {mustBeNumeric, mustBePositive, mustBeFinite} = 100000
end

arguments (Output)
    mspl struct
end

validate_spsd_struct(spsd, 'third_octave_sound_pressure_level', ...
    'required_fields', {{'fs'}}, 'check_dimensions', false);

fmax = fmax_warning(spsd.fs/2, fmax);

octave = 3;
base = 2;

mspl = band_sound_pressure_level(spsd, octave, base, fmin, fmax);
mspl.units = 'dB re 1 uPa';
mspl.name = 'Third Octave Sound Pressure Level';

end
