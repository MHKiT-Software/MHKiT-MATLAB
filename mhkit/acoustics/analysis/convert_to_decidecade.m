function out = convert_to_decidecade(spsd)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert sound pressure spectral density to decidecade spacing
%
% Parameters
% ------------
%   spsd: struct
%       Sound pressure spectral density structure containing:
%       spsd.data : Spectral density data [Pa^2/Hz or V^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector
%
% Returns
% ---------
%   out: struct
%       Decidecade band-averaged spectral density structure containing:
%       out.data : Band-averaged spectral density data [Units^2/Hz]
%       out.freq : Center frequencies of bands [Hz]
%       out.time : Time vector
%       out.name : Descriptive string
%       out.units : Units string
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
end

arguments (Output)
    out struct
end

out = convert_to_custom_bands(spsd, 10, 10, false);
out.name = 'Decidecade Sound Pressure Spectral Density';

end
