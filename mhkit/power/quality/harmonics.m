function harmonics_result = harmonics(input_data, data_sample_rate_hz, grid_freq_hz, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the harmonics from time series of voltage or current based on IEC 61000-4-7.
%
% Parameters
% ------------
%   input_data: structure
%       input_data.current : Current time series [A], n_samples x n_signals
%                            (used if present)
%       input_data.voltage : Voltage time series [V], n_samples x n_signals
%                            (used if current is not present)
%       input_data.time    : Time vector [s], n_samples x 1 (required)
%
%   data_sample_rate_hz: double
%       Sample rate of the time series [Hz]. IEC TS 62600-30 clause 7.1.4
%       requires at least 20 kHz per channel for harmonic measurements.
%
%   grid_freq_hz: double
%       Nominal grid frequency [Hz]. Options = 50 or 60
%
%   options.tolerance_percent: double (optional, default = 1)
%       Allowed deviation [%] between the sample rate implied by
%       input_data.time and data_sample_rate_hz. See mhkit_validate_sample_rate_hz.
%
% Returns
% ---------
%   harmonics_result: structure
%       harmonics_result.amplitude : Amplitude A of the sinusoid A*sin(2*pi*f*t)
%                                    at each frequency f, n_freqs x n_signals
%                                    [A] for current, [V] for voltage.
%                                    The 0 Hz row is abs(mean) of the signal.
%       harmonics_result.harmonic  : Frequencies [Hz], n_freqs x 1
%                                    = 0:5:((max_harmonic + 1) * grid_freq_hz - 5)
%       harmonics_result.type      : 'current' or 'voltage'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        input_data struct
        data_sample_rate_hz (1,1) double {mustBePositive}
        grid_freq_hz (1,1) double {mustBeMember(grid_freq_hz, [50, 60])}
        options.tolerance_percent (1,1) double {mustBePositive} = 1
    end

    % Frequency spacing, IEC 61000-4-7 clause 3.4.1 NOTE 2
    freq_spacing_hz = 5;

    % Highest harmonic order required by IEC TS 62600-30 clause 7.3
    max_harmonic = 50;

    % Validate input data structure
    if ~isfield(input_data, 'time')
        error('MHKiT:harmonics: input_data structure must contain time field');
    end

    if isfield(input_data, 'current')
        signal_data = input_data.current;
        signal_type = 'current';
    elseif isfield(input_data, 'voltage')
        signal_data = input_data.voltage;
        signal_type = 'voltage';
    else
        error('MHKiT:harmonics: input_data structure must contain either current or voltage field');
    end

    if ~isnumeric(signal_data) || ~isreal(signal_data)
        error('MHKiT:harmonics: %s data must be a real numeric array', signal_type);
    end

    % Validate time vector sample rate
    sample_rate_validation = mhkit_validate_sample_rate_hz(input_data.time, data_sample_rate_hz, ...
        'tolerance_percent', options.tolerance_percent);

    if ~sample_rate_validation.pass
        error(['MHKiT:harmonics: Time vector sample rate validation failed.\n' ...
               'Expected: %.2f Hz, Observed: %.2f Hz (median), Deviation: %.2f%%, Tolerance: %.2f%%\n' ...
               'Time format detected: %s'], ...
               data_sample_rate_hz, sample_rate_validation.median_sample_rate_hz, ...
               sample_rate_validation.deviation_percent, sample_rate_validation.tolerance_percent, ...
               sample_rate_validation.time_format);
    end

    % Validate dimensions of every signal field against the time vector
    num_time_samples = length(input_data.time);
    field_names = fieldnames(input_data);

    for field_idx = 1:length(field_names)
        field_name = field_names{field_idx};
        if strcmp(field_name, 'time')
            continue;
        end
        field_rows = size(input_data.(field_name), 1);
        if field_rows ~= num_time_samples
            error('MHKiT:harmonics: %s data rows (%d) must match time vector length (%d)', ...
                field_name, field_rows, num_time_samples);
        end
    end

    num_signal_columns = size(signal_data, 2);

    % Define frequencies [Hz] for harmonic amplitude reporting: 0, 5, 10, ... up to
    % just below harmonic order max_harmonic + 1 (3055 Hz for a 60 Hz grid,
    % 2545 Hz for a 50 Hz grid).
    %
    % 5 Hz spacing matches IEC 61000-4-7, which analyzes 200 ms windows
    % (1 / 0.2 s = 5 Hz).
    %
    % Frequencies above harmonic order max_harmonic are included because
    % harmonic_subgroups and interharmonics use them for order max_harmonic.
    max_freq_hz = (max_harmonic + 1) * grid_freq_hz - freq_spacing_hz;
    harmonic_freq_grid_hz = (0:freq_spacing_hz:max_freq_hz)';
    num_freqs = length(harmonic_freq_grid_hz);

    % FFT amplitude of every column
    fft_amplitude = abs(fft(signal_data, [], 1));

    % Map each output frequency to the nearest FFT bin: bin k has frequency k * fs / n
    exact_bin_index = harmonic_freq_grid_hz * num_time_samples / data_sample_rate_hz;
    closest_bin_idx = round(exact_bin_index) + 1;  % +1 for MATLAB 1-based indexing
    nyquist_bin_idx = floor(num_time_samples / 2) + 1;
    nyquist_freq_hz = data_sample_rate_hz / 2;
    in_band = harmonic_freq_grid_hz <= nyquist_freq_hz & closest_bin_idx <= nyquist_bin_idx;

    if ~all(in_band)
        warning('MHKiT:harmonics:AboveNyquist', ...
            ['Frequencies from %.1f Hz upward are above the Nyquist frequency of the ' ...
             '%.1f Hz sample rate and are set to zero. IEC TS 62600-30 clause 7.1.4 ' ...
             'requires at least 20 kHz per channel.'], ...
            min(harmonic_freq_grid_hz(~in_band)), data_sample_rate_hz);
    end

    harmonics_reindexed = zeros(num_freqs, num_signal_columns);
    harmonics_reindexed(in_band, :) = fft_amplitude(closest_bin_idx(in_band), :);

    % Single-sided normalization: 2/n for all frequencies, 1/n for DC
    harmonics_normalized = harmonics_reindexed / num_time_samples * 2;
    harmonics_normalized(1, :) = harmonics_reindexed(1, :) / num_time_samples;

    % Create output structure
    harmonics_result = struct();
    harmonics_result.amplitude = harmonics_normalized;
    harmonics_result.harmonic = harmonic_freq_grid_hz;
    harmonics_result.type = signal_type;

end
