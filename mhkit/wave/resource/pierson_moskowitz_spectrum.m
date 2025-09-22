function S = pierson_moskowitz_spectrum(frequency, Tp, Hs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Pierson-Moskowitz Spectrum from IEC TS 62600-2 ED2 Annex C.2 (2019)
%
% Parameters
% ------------
%     frequency: vector
%         Wave frequency [Hz]
%
%     Tp: float
%         Peak period [s]
%
%     Hs: float
%         Significant wave height [m]
%
% Returns
% ---------
%     S: structure with fields
%         .spectrum  = Spectral density [m^2/Hz]
%         .frequency = Frequency [Hz]
%         .type      = 'Pierson-Moskowitz (Hs,Tp)'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    frequency {mustBeNumeric, mustBeFinite}
    Tp {mustBeNumeric, mustBeFinite, mustBePositive}
    Hs {mustBeNumeric, mustBeFinite, mustBePositive}
end

% Use JONSWAP spectrum with gamma = 1 (equivalent to Pierson-Moskowitz)
S_jonswap = jonswap_spectrum(frequency, Tp, Hs, 1);

% Create output structure with Pierson-Moskowitz naming for backward compatibility
S = struct();
S.frequency = S_jonswap.frequency;
S.spectrum = S_jonswap.spectrum;
S.type = sprintf('Pierson-Moskowitz (%.1fm,%.1fs)', Hs, Tp);
