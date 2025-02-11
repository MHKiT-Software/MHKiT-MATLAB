function S=jonswap_spectrum(frequency,Tp,Hs,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates JONSWAP spectrum from Hasselmann et al (1973)
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
%         Significant Wave Height (s)
%
%     gamma: float (optional)
%         Gamma
%
% Returns
% ---------
%     S: structure
%
%
%         S.spectrum=Spectral Density (m^2/Hz)
%
%         S.type=String of the spectra type, i.e. Bretschneider,
%         time series, date stamp etc.
%
%         S.frequency= frequency (Hz)
%
%         S.Te
%
%         S.Hm0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isa(frequency,'py.numpy.ndarray') ~= 1)
    frequency = py.numpy.array(frequency);
end

if nargin == 3
        S_py=py.mhkit.wave.resource.jonswap_spectrum(frequency,Tp,Hs);
elseif nargin == 4
        S_py=py.mhkit.wave.resource.jonswap_spectrum(frequency, Tp, Hs, pyargs('gamma', varargin{1}));
end

S_py = typecast_from_mhkit_python(S_py);

S = struct();

S.frequency = S_py.index.data;
S.spectrum = S_py.data;
S.type = S_py.columns{1};
