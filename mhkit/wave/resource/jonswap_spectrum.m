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

if isrow(frequency)
    f = sort(frequency');
else
    f = sort(frequency);
end

B_PM = (5/4)*(1/Tp)^4;
A_PM = B_PM*(Hs/2)^2;
S_f  = A_PM.*f.^(-5).*exp(-B_PM.*f.^(-4));

if nargin == 3
    TpsqrtHs = Tp/np.sqrt(Hs);
    if TpsqrtHs <= 3.6
            gamma = 5;
    elseif TpsqrtHs > 5
        gamma = 1;
    else
        gamma = np.exp(5.75 - 1.15*TpsqrtHs);
    end
else
    gamma = varargin{1};
end

% Cutoff frequencies for gamma function
siga = 0.07;
sigb = 0.09;

fp = 1/Tp; % peak frequency
lind = f<=fp;
hind = f>fp;
Gf = zeros(size(f));
Gf(lind) = gamma.^exp(-(f(lind)-fp).^2/(2*siga^2*fp.^2));
Gf(hind) = gamma.^exp(-(f(hind)-fp).^2/(2*sigb^2*fp.^2));
C = 1- 0.287*log(gamma);
Sf = C.*S_f.*Gf;

S.spectrum = Sf;
S.frequency = f;
S.Tp = Tp;
