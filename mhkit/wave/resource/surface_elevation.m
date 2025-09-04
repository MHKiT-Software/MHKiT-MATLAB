function wave_elevation = surface_elevation(S, time_index, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate wave surface elevation from spectrum
%
% Parameters
% ------------
%     S: structure
%         Spectral data with fields:
%             S.spectrum: (n_freq x 1) spectral density values (m^2/Hz)
%             S.frequency: (n_freq x 1) frequency vector (Hz)
%     
%     time_index: vector
%         (n_time x 1) time vector (s)
%
%     options: structure (optional)
%         Optional fields:
%             seed: random seed (default = 123)
%             frequency_bins: bin widths for each frequency
%             phases: explicit phases (radians)
%             method: 'ifft' (default) or 'sum_of_sines'
%
% Returns
% ---------
%     wave_elevation: structure
%         Generated wave elevation with fields:
%             elevation: (n_time x 1) wave elevation (m)
%             time: time vector (s)
%             type: description string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S struct
    time_index (:,1) double
    options.seed (1,1) double = 123
    options.frequency_bins = []
    options.phases = []
    options.method char {mustBeMember(options.method, {'ifft','sum_of_sines'})} = 'ifft'
end

% Extract frequency and spectrum
f = S.frequency(:);
Sf = S.spectrum(:);
Nf = numel(f);
Nt = numel(time_index);

% Handle frequency bins
if isempty(options.frequency_bins)
    df = diff(f);
    df_uniform = all(abs(df - df(1)) < 1e-8);
    df_val = mean(df);
    delta_f = df_val * ones(size(f));
else
    delta_f = options.frequency_bins(:);
    if length(delta_f) ~= length(f)
        error('frequency_bins must match the length of frequency vector.');
    end
    df_val = delta_f(1);
    df_uniform = all(abs(delta_f - df_val) < 1e-8);
end

% Handle phases
if isempty(options.phases)
    rng(options.seed);
    phase = 2*pi*rand(Nf, 1);
else
    phase = options.phases(:);
    if length(phase) ~= Nf
        error('phases must match the length of frequency vector.');
    end
end

% Choose method
omega = 2*pi*f;
A = sqrt(2 * Sf .* delta_f);

switch options.method
    case 'ifft'
        if f(1) ~= 0 || ~df_uniform
            warning('Switching to sum_of_sines because ifft requires f(1)==0 and uniform spacing.');
            options.method = 'sum_of_sines';
        else
            % Create complex spectrum
            A_complex = A .* exp(1i * phase);
            % Use irfft approximation: 0 to fN mapped to 0 to Nyquist
            % Double spectrum (except DC and Nyquist if exists)
            spectrum_full = [A_complex(1); A_complex(2:end-1)/2; conj(flip(A_complex(2:end-1)/2))];
            eta = real(ifft(spectrum_full, Nt)) * Nt;
            wave_elevation.elevation = eta;
    end

    case 'sum_of_sines'
        B = omega .* time_index';
        B = B'; % (Nt x Nf)
        C = cos(B + phase');
        eta = C * A;
        wave_elevation.elevation = eta;
end

wave_elevation.time = time_index;
wave_elevation.type = 'Time Series from Spectra';
end
