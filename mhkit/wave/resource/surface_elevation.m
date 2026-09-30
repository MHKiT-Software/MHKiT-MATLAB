function wave_elevation = surface_elevation(S, time_index, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave elevation time-series from spectrum
%
% Parameters
% ------------
% S : struct
%   Wave spectrum structure:
%     S.spectrum : vector [m^2/Hz]
%       Spectral density
%     S.frequency : vector [Hz]
%       Frequency
% time_index : vector [s]
%   Time used to create the wave elevation time-series,
%   for example, time_index = 0:0.01:100
% seed : double (optional)
%   Random seed. Default = [] (unseeded, non-reproducible)
% frequency_bins : vector [Hz] (optional)
%   Bin widths for frequency of S. Required for unevenly sized bins.
% phases : vector [rad] (optional)
%   Explicit phases for frequency components (overrides seed),
%   for example, phases = rand(length(S.frequency), 1) * 2 * pi
% method : char (optional)
%   Method used to calculate the surface elevation. 'ifft' (Inverse
%   Fast Fourier Transform) used by default if the given frequency_bins
%   is empty or evenly spaced. 'sum_of_sines' explicitly sums each
%   frequency component and is used by default if uneven frequency_bins
%   are provided. The 'ifft' method is significantly faster.
%   Default = 'ifft'
%
% Returns
% ---------
% wave_elevation : struct
%   wave_elevation.elevation : vector [m]
%     Wave surface elevation
%   wave_elevation.time : vector [s]
%     Time vector
%   wave_elevation.type : char
%     'Time Series from Spectra'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    S struct
    time_index (:,1) double
    options.seed {mustBeNumeric} = []
    options.frequency_bins = []
    options.phases = []
    options.method char {mustBeMember(options.method, {'ifft','sum_of_sines'})} = 'ifft'
end

arguments (Output)
    wave_elevation struct
end

if ~isempty(options.seed) && ~isscalar(options.seed)
    error('MHKiT:surface_elevation:InvalidInput', 'seed must be a scalar or empty (unseeded).');
end

% Extract frequency and spectrum
f = S.frequency(:);
Sf = S.spectrum(:);
Nf = numel(f);
Nt = numel(time_index);

if isempty(options.frequency_bins)
    delta_f = f(2) - f(1);
    df = diff(f);
    df_uniform = all(abs(df - df(1)) < 1e-8);
else
    freq_bins = options.frequency_bins(:);
    if length(freq_bins) ~= length(f)
        error('frequency_bins must match the length of frequency vector.');
    end
    df_uniform = all(abs(freq_bins - freq_bins(1)) < 1e-8);
    if df_uniform
        delta_f = freq_bins(1);
    else
        delta_f = freq_bins;
    end
end

% An empty seed is left unseeded (non-reproducible),
% following MHKiT-Python convention.
% https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L372-L373
if isempty(options.phases)
    if ~isempty(options.seed)
        rng(options.seed);
    end
    phase = 2*pi*rand(Nf, 1);
else
    phase = options.phases(:);
    if length(phase) ~= Nf
        error('phases must match the length of frequency vector.');
    end
end

method = options.method;
if strcmp(method, 'ifft')
    if f(1) ~= 0
        warning('MHKiT:surface_elevation:MethodFallback', ...
            ['ifft method must have zero frequency defined. Setting ' ...
            'method to less efficient sum_of_sines method.']);
        method = 'sum_of_sines';
    end
    if ~df_uniform
        warning('MHKiT:surface_elevation:MethodFallback', ...
            ['ifft method must have evenly spaced frequency bins. ' ...
            'Setting method to less efficient sum_of_sines method.']);
        method = 'sum_of_sines';
    end
end

omega = 2*pi*f;
A = sqrt(2 * Sf .* delta_f);

if strcmp(method, 'ifft')
    A_complex = A .* (cos(phase) + 1i * sin(phase));
    A_scaled = 0.5 * A_complex * length(time_index);
    % Use MATLAB's equivalent of irfft
    eta = real(ifft(A_scaled, length(time_index), 'symmetric'));
    wave_elevation.elevation = eta;
elseif strcmp(method, 'sum_of_sines')
    B = omega .* time_index';
    B = B'; % (Nt x Nf)
    C = cos(B + phase');
    eta = C * A;
    wave_elevation.elevation = eta;
end

wave_elevation.time = time_index;
wave_elevation.type = 'Time Series from Spectra';
end
