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

% DEBUG: Print all dimensions before calculation
fprintf('DEBUG - Variable dimensions in mler_coefficients:\n');
fprintf('  RAO: %s\n', mat2str(size(RAO)));
fprintf('  wave_spectrum: %s\n', mat2str(size(wave_spectrum)));
fprintf('  freq: %s\n', mat2str(size(freq)));
fprintf('  dw: %s (scalar: %g)\n', mat2str(size(dw)), dw);
fprintf('  m0: %s (value: %g)\n', mat2str(size(m0)), m0);
fprintf('  m1: %s (value: %g)\n', mat2str(size(m1)), m1);
fprintf('  m2: %s (value: %g)\n', mat2str(size(m2)), m2);
fprintf('  wBar: %s (value: %g)\n', mat2str(size(wBar)), wBar);

% calculate coefficient_a from Quon2016 Eqn.8
% Explicit scalar handling for MATLAB R2022b/R2023b compatibility
% Debug output before computation
fprintf('  BEFORE explicit scalar conversion:\n');
fprintf('    m0 original: size=%s, class=%s, value=%g\n', mat2str(size(m0)), class(m0), m0);
fprintf('    m1 original: size=%s, class=%s, value=%g\n', mat2str(size(m1)), class(m1), m1);
fprintf('    m2 original: size=%s, class=%s, value=%g\n', mat2str(size(m2)), class(m2), m2);
fprintf('    wBar original: size=%s, class=%s, value=%g\n', mat2str(size(wBar)), class(wBar), wBar);

% Ensure all scalars are explicitly converted to same size as vectors
m0_scalar = double(m0(1));  % Ensure scalar
m1_scalar = double(m1(1));  % Ensure scalar
m2_scalar = double(m2(1));  % Ensure scalar
wBar_scalar = double(wBar(1));  % Ensure scalar

fprintf('  AFTER explicit scalar conversion:\n');
fprintf('    m0_scalar: size=%s, class=%s, value=%g\n', mat2str(size(m0_scalar)), class(m0_scalar), m0_scalar);
fprintf('    m1_scalar: size=%s, class=%s, value=%g\n', mat2str(size(m1_scalar)), class(m1_scalar), m1_scalar);
fprintf('    m2_scalar: size=%s, class=%s, value=%g\n', mat2str(size(m2_scalar)), class(m2_scalar), m2_scalar);
fprintf('    wBar_scalar: size=%s, class=%s, value=%g\n', mat2str(size(wBar_scalar)), class(wBar_scalar), wBar_scalar);

% Calculate each term explicitly with size monitoring
term1 = abs(RAO);  % [1 x N]
fprintf('  term1 = abs(RAO): %s\n', mat2str(size(term1)));

term2 = sqrt(2 .* dw .* wave_spectrum);  % [1 x N]
fprintf('  term2 = sqrt(2.*dw.*wave_spectrum): %s\n', mat2str(size(term2)));

term3a = m2_scalar - (freq .* m1_scalar);  % [1 x N]
fprintf('  term3a = m2_scalar - (freq.*m1_scalar): %s\n', mat2str(size(term3a)));

term3b = (freq .* m0_scalar) - m1_scalar;  % [1 x N]
fprintf('  term3b = (freq.*m0_scalar) - m1_scalar: %s\n', mat2str(size(term3b)));

term3c = wBar_scalar .* term3b;  % [1 x N]
fprintf('  term3c = wBar_scalar .* term3b: %s\n', mat2str(size(term3c)));

term3 = term3a + term3c;  % [1 x N]
fprintf('  term3 = term3a + term3c: %s\n', mat2str(size(term3)));

% Final calculation with explicit element-wise operations
numerator = term1 .* term2 .* term3;  % [1 x N]
fprintf('  numerator = term1 .* term2 .* term3: %s\n', mat2str(size(numerator)));

denominator_scalar = (m0_scalar * m2_scalar) - (m1_scalar ^ 2);  % scalar
fprintf('  denominator_scalar: size=%s, value=%g\n', mat2str(size(denominator_scalar)), denominator_scalar);

coeff_a_rn = numerator ./ denominator_scalar;  % [1 x N]
fprintf('  coeff_a_rn final: %s\n', mat2str(size(coeff_a_rn)));

% phase delay should be positive number
phase = unwrap(angle(RAO));
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
