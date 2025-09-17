function spsdl = sound_pressure_spectral_density_level(spsd)
% 
% Calculates the sound pressure spectral density level from
% the mean square sound pressure spectral density.
% 
% Parameters
% ----------
% spsd: struct
%     Mean square sound pressure spectral density in Pa^2/Hz
% 
% Returns
% -------
% spsdl: struct
%     Sound pressure spectral density level [dB re 1 uPa^2/Hz]
%     indexed by time and frequency

arguments (Input)
    spsd struct
end

arguments (Output)
    spsdl struct
end

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