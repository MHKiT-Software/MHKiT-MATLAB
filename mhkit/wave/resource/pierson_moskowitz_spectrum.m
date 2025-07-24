function S = pierson_moskowitz_spectrum(frequency, Tp, Hs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates Pierson-Moskowitz Spectrum from Tucker and Pitt (2001)
%
% Parameters
% ------------
%
%     Frequency: float
%         Wave frequency (Hz)
%
%     Tp: float
%         Peak Period (s)
%
%     Hs: float
%         Significant wave height (m)
%
%
% Returns
% ---------
%     S: structure
%
%         S.spectrum=Spectral Density (m^2/Hz)
%
%         S.type=String of the spectra type, i.e. (Pierson-Moskowitz 8.0s)
%
%         S.frequency= frequency (Hz)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isa(frequency,'py.numpy.ndarray') ~= 1)
    frequency = py.numpy.array(frequency);
end

S_py = py.mhkit.wave.resource.pierson_moskowitz_spectrum(frequency, Tp, Hs);

S_py = typecast_from_mhkit_python(S_py);

S = struct();

S.frequency = S_py.index.data;
S.spectrum = S_py.data;
S.type = S_py.columns{1};
