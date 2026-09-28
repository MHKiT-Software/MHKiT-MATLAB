function mler = mler_coefficients(RAO, wave_spectrum, response_desired)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculate MLER (most likely extreme response) coefficients from a
%     sea state spectrum and a response RAO.
%
%     Parameters
%     ----------
%         RAO : array (N x 1)
%             Response amplitude operator [-]
%         wave_spectrum: struct
%             wave_spectrum.spectrum - Spectral density [m^2/Hz] (N x 1)
%             wave_spectrum.frequency - Frequency [Hz] (N x 1)
%         response_desired: scalar
%             Desired response amplitude
%
%      Returns
%      -------
%         mler : struct
%             mler.conditioned_spectrum - Conditioned wave spectral amplitude [m^2-s] (N x 1)
%             mler.phase - Phase [rad] (N x 1)
%             mler.frequency - Frequency [Hz] (N x 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(RAO,'numeric')
    error('ERROR: rao must be a double array')
end
if ~isa(wave_spectrum,'struct')
    error('ERROR: wave_spectrum must be struct with fieldnames spectrum and frequency')
end
if ~isa(response_desired,'numeric')
    error('ERROR: response_desired must be an int or double')
end

% validate all inputs have same length
N = length(RAO);
if length(wave_spectrum.frequency) ~= N || length(wave_spectrum.spectrum) ~= N
    error('MHKiT:loads:mler_coefficients: RAO, frequency, and spectrum must have same length');
end

% convert from Hz to rad/s
freq_rad = wave_spectrum.frequency * (2*pi);
wave_spectrum_rad = wave_spectrum.spectrum / (2*pi);
dw = (2*pi - 0) / (N-1);

% response spectrum
R.spectrum = abs(RAO).^2 .* (2*wave_spectrum_rad);
R.type = 'response';
R.frequency = freq_rad;

% spectral moment calculations
m0 = frequency_moment(R, 0);
m1 = frequency_moment(R, 1);
m2 = frequency_moment(R, 2);
wBar = m1/m0;

% calculate coefficient_a from Quon2016 Eqn.8
coeff_a_rn = abs(RAO) .* sqrt(2*dw.*wave_spectrum_rad) .* ((m2 - freq_rad.*m1) + wBar.*(freq_rad.*m0 - m1)) ./ (m0*m2 - m1^2);
% phase delay should be positive number
phase = -unwrap(angle(RAO));
% for negative values of Amp, add pi phase shift, flip sign
phase(coeff_a_rn < 0) = phase(coeff_a_rn < 0) - pi;
coeff_a_rn(coeff_a_rn < 0) = coeff_a_rn(coeff_a_rn < 0) * -1;

% calculate conditioned spectrum [m^2-s/rad]
conditioned_spectrum = wave_spectrum_rad .* coeff_a_rn.^2 .* response_desired^2;
conditioned_spectrum(isnan(conditioned_spectrum)) = 0; % replace nans with zero
% if the response amplitude we ask for is negative, we will add
% a pi phase shift to the phase information.  This is because
% the sign of response_desired is lost in the squaring above.
% Ordinarily this would be put into the final equation, but we
% are shaping the wave information so that it is buried in the
% new spectral information, S. (AP)
if response_desired < 0
    phase = phase + pi;
end

% outputs
mler.conditioned_spectrum = conditioned_spectrum;
mler.phase = phase;
mler.frequency = wave_spectrum.frequency;

end
