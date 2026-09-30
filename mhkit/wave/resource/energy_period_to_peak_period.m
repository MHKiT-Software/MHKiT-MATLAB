function Tp = energy_period_to_peak_period(Te, gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert spectral energy period (Te) to peak period (Tp)
%
% Uses ITTC approximation for JONSWAP Spectrum from "The Specialist
% Committee on Waves, Final Report and Recommendations to the 23rd ITTC",
% Proceedings of the 23rd ITTC - Volume 2, Table A4.
%
% Parameters
% ------------
% Te : numeric [s]
%   Spectral energy period (scalar, vector, or array)
% gamma : double
%   Peak enhancement factor for JONSWAP spectrum
%
% Returns
% ---------
% Tp : numeric [s]
%   Spectral peak period (same shape as Te)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Te {mustBeNumeric}
    gamma (1,1) {mustBeNumeric, mustBePositive}
end

arguments (Output)
    Tp {mustBeNumeric}
end

% ITTC approximation factor
factor = 0.8255 + 0.03852 * gamma - 0.005537 * gamma^2 + 0.0003154 * gamma^3;

Tp = Te / factor;

end
