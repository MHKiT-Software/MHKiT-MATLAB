function spsdl = sound_pressure_spectral_density_level(spsd)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound pressure spectral density level from mean square sound pressure spectral density
%
% Parameters
% ------------
%   spsd: struct
%       spsd.data : Mean square sound pressure spectral density [Pa^2/Hz]
%       spsd.freq : Frequency vector [Hz] (optional)
%       spsd.time : Time vector (optional)
%
% Returns
% ---------
%   spsdl: struct
%       spsdl.data : Sound pressure spectral density level [dB re 1 uPa^2/Hz]
%       spsdl.freq : Frequency vector [Hz] (if present in input)
%       spsdl.time : Time vector (if present in input)
%       spsdl.name : Data name identifier
%       spsdl.units : Data units string
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
end

arguments (Output)
    spsdl struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'sound_pressure_spectral_density_level', ...
    'require_positive', true, 'check_dimensions', false);

% reference value of sound pressure
reference = 1e-12; % Pa^2 to 1uPa^2

% Sound pressure spectral density level from mean square values
lpf = 10 * log10(spsd.data / reference);

% update struct
spsdl = spsd;
spsdl.data = lpf;
spsdl.name = "Sound Pressure Spectral Density Level";
spsdl.units = "dB re 1 uPa^2/Hz";

end
