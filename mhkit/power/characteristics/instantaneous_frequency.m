function frequency = instantaneous_frequency(voltage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate instantaneous frequency of measured voltage
%
% Parameters
% ------------
% voltage : struct
%   Measured voltage time series
%     voltage.voltage : double [V]
%       Measured voltage, n_samples x n_signals, one time series per column
%     voltage.time : double [s]
%       Time vector, n_samples x 1
%
% Returns
% ---------
% frequency : struct
%   Instantaneous frequency of each voltage signal
%     frequency.frequency : double [Hz]
%       Instantaneous frequency, (n_samples - 1) x n_signals
%     frequency.time : double [s]
%       Time vector, (n_samples - 1) x 1. One element shorter than the
%       input because the phase is differentiated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        voltage struct
    end

    % Validate input structure has required fields
    if ~isfield(voltage, 'voltage')
        error('MHKiT:instantaneous_frequency:InvalidInput', 'voltage structure must contain voltage field');
    end
    if ~isfield(voltage, 'time')
        error('MHKiT:instantaneous_frequency:InvalidInput', 'voltage structure must contain time field');
    end

    % Extract data from structure
    voltage_data = voltage.voltage;
    time_vector = voltage.time;

    % Validate dimensions
    if size(voltage_data, 1) ~= length(time_vector)
        error('MHKiT:instantaneous_frequency:InvalidInput', 'voltage data rows must match time vector length');
    end

    % Get data dimensions
    [num_samples, num_columns] = size(voltage_data);

    % Validate minimum data length for meaningful frequency calculation
    if num_samples < 4
        error('MHKiT:instantaneous_frequency:InvalidInput', 'voltage data must have at least 4 samples for frequency calculation');
    end

    % Warn the user if the sample interval varies by more than the tolerance.
    % The phase derivative below divides by each local dt, but the FFT-based
    % Hilbert transform that produces the phase assumes uniform sampling.
    % Irregular intervals corrupt the phase itself before dt is applied, so
    % using the local dt cannot correct for them.
    sample_rate = mhkit_validate_sample_rate_hz(time_vector);
    if ~sample_rate.is_uniform
        warning('MHKiT:instantaneous_frequency:SampleRateVariation', ...
                ['The sample interval of this signal varies by more than %g%% from the mean. ', ...
                 'The FFT-based Hilbert transform assumes uniform sampling, so the ', ...
                 'instantaneous phase, and therefore the instantaneous frequency, is ', ...
                 'likely to be inaccurate. ', ...
                 'Mean sample rate: %g Hz, max: %g Hz, min: %g Hz, standard deviation: %g Hz'], ...
                sample_rate.tolerance_percent, sample_rate.mean_sample_rate_hz, ...
                sample_rate.max_sample_rate_hz, sample_rate.min_sample_rate_hz, ...
                sample_rate.std_sample_rate_hz);
    end

    % Calculate time differences for frequency calculation
    time_diff = diff(time_vector(:));  % Ensure column vector

    % Analytic signal of every column
    analytic_signal = mhkit_hilbert(voltage_data);

    % Instantaneous phase with 2*pi discontinuities removed
    unwrapped_phase = unwrap(angle(analytic_signal), [], 1);

    % Instantaneous frequency
    frequency_data = diff(unwrapped_phase, 1, 1) ./ (2.0 * pi * time_diff);

    % Create output structure
    frequency = struct();
    frequency.frequency = frequency_data;
    frequency.time = time_vector(2:end);

end

