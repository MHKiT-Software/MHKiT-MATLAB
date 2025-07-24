function mler = mler_coefficients(RAO, wave_spectrum, response_desired)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculate MLER (most likely extreme response) coefficients from a
%     sea state spectrum and a response RAO.
%
%     Parameters
%     ----------
%         RAO : array
%             Response amplitude operator [-]
%         wave_spectrum: struct
%             Struct with wave spectral density [m^2/Hz] and frequency [Hz]
%         response_desired: int or float
%             Latitude longitude pairs at which to extract data.
%
%      Returns
%      -------
%         mler : struct
%             containing conditioned wave spectral amplitude
%             coefficient [m^2-s], and phase [rad] indexed by frequency [Hz].
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

% convert from Hz to rad/s
freq = wave_spectrum.frequency * (2*pi);
freq_hz = wave_spectrum.frequency;
wave_spectrum = wave_spectrum.spectrum / (2*pi);
dw = (2*pi - 0) / (length(freq)-1);

% response spectrum
R.spectrum = abs(RAO).^2 .* (2*wave_spectrum);
R.type = 'response';
R.frequency = freq;

% spectral moment calculations
m0 = frequency_moment(R, 0);
m1 = frequency_moment(R, 1);
m2 = frequency_moment(R, 2);
wBar = m1/m0;

% calculate coefficient_a from Quon2016 Eqn.8
coeff_a_rn = abs(RAO) .* sqrt(2*wave_spectrum*dw) .* ((m2 - freq*m1) + wBar*(freq*m0 - m1)) ./ (m0*m2 - m1^2);
% phase delay should be positive number
phase = -unwrap(angle(RAO));
% for negative values of Amp, add pi phase shift, flip sign
phase(coeff_a_rn < 0) = phase(coeff_a_rn < 0) - pi;
coeff_a_rn(coeff_a_rn < 0) = coeff_a_rn(coeff_a_rn < 0) * -1;

% calculate conditioned spectrum [m^2-s/rad]
S = wave_spectrum .* coeff_a_rn.^2 .* response_desired^2;
S(isnan(S)) = 0; % replace nans with zero
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
mler.conditioned_spectrum = S;
mler.phase = phase;
mler.frequency = freq_hz;

end
