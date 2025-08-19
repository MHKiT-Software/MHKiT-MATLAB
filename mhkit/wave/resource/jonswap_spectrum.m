function S = jonswap_spectrum(frequency, Tp, Hs, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates JONSWAP spectrum based on IEC TS 62600-2 ED2 Annex C.2 (2019)
%
% Parameters
% ------------
%     frequency: vector
%         Wave frequency (Hz)
%
%     Tp: float
%         Peak Period (s)
%
%     Hs: float
%         Significant Wave Height (m)
%
%     gamma: float (optional)
%         Peak enhancement factor
%
% Returns
% ---------
%     S: structure with fields
%         .spectrum  = Spectral Density (m^2/Hz)
%         .frequency = Frequency (Hz)
%         .type      = 'JONSWAP (Hs,Tp)'
%         .Hm0       = Significant wave height (m)
%         .Te        = Energy period (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure column vector
frequency = frequency(:);
f = frequency;
df = diff(f);
uniform_spacing = all(abs(df - df(1)) < 1e-6);
if uniform_spacing
    df_val = df(1);
else
    df_val = [diff(f); diff(f(end-1:end))]; % Extend last bin
end

% Constants
fp = 1 / Tp;
B_PM = (5 / 4) * (1 / Tp)^4;
A_PM = B_PM * (Hs / 2)^2;

% Initialize spectral density
S_f = zeros(size(f));
nonzero = f > 0;
S_f(nonzero) = A_PM .* f(nonzero).^(-5) .* exp(-B_PM .* f(nonzero).^(-4));

% Gamma computation
if nargin < 4 || isempty(gamma)
    TpsqrtHs = Tp / sqrt(Hs);
    if TpsqrtHs <= 3.6
        gamma = 5;
    elseif TpsqrtHs > 5
        gamma = 1;
    else
        gamma = exp(5.75 - 1.15 * TpsqrtHs);
    end
end

% Spreading function G(f)
siga = 0.07;
sigb = 0.09;
Gf = zeros(size(f));
lind = f <= fp;
hind = f > fp;
Gf(lind) = gamma .^ exp(-((f(lind) - fp).^2) ./ (2 * siga^2 * fp^2));
Gf(hind) = gamma .^ exp(-((f(hind) - fp).^2) ./ (2 * sigb^2 * fp^2));

C = 1 - 0.287 * log(gamma);
Sf = C * S_f .* Gf;

% Outputs
S.frequency = f;
S.spectrum = Sf;
S.type = sprintf('JONSWAP (%.2fm, %.2fs)', Hs, Tp);

% Compute Hm0 and Te
if uniform_spacing
    m0 = sum(Sf) * df_val;
    m1 = sum(Sf ./ f) * df_val;
else
    m0 = sum(Sf .* df_val);
    m1 = sum((Sf ./ f) .* df_val);
end

S.Hm0 = 4 * sqrt(m0);
S.Te = m1 / m0;

end
