function S = jonswap_spectrum(frequency, Tp, Hs, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates JONSWAP spectrum from IEC TS 62600-2 ED2 Annex C.2 (2019)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure column vector and sort (match MHKiT-Python behavior)
frequency = frequency(:);
f = sort(frequency);

% Constants (same as MHKiT-Python)
B_PM = (5 / 4) * (1 / Tp)^4;
A_PM = B_PM * (Hs / 2)^2;

% Avoid divide by zero if the 0 frequency is provided
% The zero frequency should always have 0 amplitude, otherwise
% we end up with a mean offset when computing the surface elevation.
S_f = zeros(size(f));
if f(1) == 0.0
    inds = 2:length(f);  % Skip first element like MHKiT-Python
else
    inds = 1:length(f);  % Use all elements
end

S_f(inds) = A_PM * f(inds).^(-5) .* exp(-B_PM * f(inds).^(-4));

% Gamma computation (match MHKiT-Python logic)
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

% Cutoff frequencies for gamma function (match MHKiT-Python)
siga = 0.07;
sigb = 0.09;

% Peak frequency
fp = 1 / Tp;
lind = f <= fp;
hind = f > fp;
Gf = zeros(size(f));
Gf(lind) = gamma .^ exp(-((f(lind) - fp).^2) ./ (2 * siga^2 * fp^2));
Gf(hind) = gamma .^ exp(-((f(hind) - fp).^2) ./ (2 * sigb^2 * fp^2));

C = 1 - 0.287 * log(gamma);
Sf = C * S_f .* Gf;

% Output structure (match MHKiT-Python naming style)
name = sprintf('JONSWAP (%.1fm,%.1fs)', Hs, Tp);
S.frequency = f;
S.spectrum = Sf;
S.type = name;

end
