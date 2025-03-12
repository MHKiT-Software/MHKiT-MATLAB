function frequency = instantaneous_frequency(measured_voltage, time_dimension)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculates instantaneous frequency of measured voltage signals
%     using Hilbert transform to extract phase information. Works with
%     both single and multi-channel voltage measurements.
%
% Parameters
% ------------
%   measured_voltage: structure
%       measured_voltage.voltage : Matrix of voltage measurements [V]
%                                 Each row is a time point, each column is a channel
%       measured_voltage.time : Time vector [s]
%   time_dimension: string (optional)
%       Currently not used in MATLAB implementation
%       Reserved for future compatibility
%
% Returns
% ---------
%   frequency: structure
%       frequency.frequency : Matrix of calculated instantaneous frequencies [Hz]
%                            Has one fewer row than the input voltage
%       frequency.time : Time vector corresponding to frequency values
%
% Key Equations
% -------------
% 1. Analytic Signal via Hilbert Transform:
%    z(t) = x(t) + j·H[x(t)]
%    where H[x(t)] is the Hilbert transform of x(t)
%
% 2. Instantaneous Phase:
%    φ(t) = unwrap(angle(z(t)))
%
% 3. Instantaneous Frequency:
%    f(t) = (1/2π) · dφ(t)/dt
%
% Notes
% ------
% The function implements a custom Hilbert transform without requiring
% the Signal Processing Toolbox.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input validation
    if ~isstruct(measured_voltage)
        error('measured_voltage must be a structure with voltage and time fields');
    end

    % Validate required fields exist
    required_fields = {'voltage', 'time'};
    if ~all(isfield(measured_voltage, required_fields))
        error('measured_voltage structure must contain fields: voltage and time');
    end

    % Extract voltage and time data
    voltage_data = measured_voltage.voltage;
    time = measured_voltage.time;

    % Validate data types
    if ~isnumeric(voltage_data) || ~isnumeric(time)
        error('Voltage and time data must be numeric');
    end

    % Check for empty data
    if isempty(voltage_data) || isempty(time)
        error('Voltage and time data cannot be empty');
    end

    % Verify time vector is suitable for frequency calculation
    if length(time) <= 1
        error('Time vector must contain at least two points for frequency calculation');
    end

    % Calculate time step
    dt = diff(time);
    if std(dt)/mean(dt) > 0.01  % Check if sampling is approximately uniform
        warning('Time steps vary by more than 1%. Frequency calculation assumes uniform sampling.');
    end
    dt = dt(1); % Assuming uniform sampling

    % Initialize output structure
    frequency = struct();

    % Get size of voltage data
    [num_samples, num_channels] = size(voltage_data);
    freq_data = zeros(num_samples-1, num_channels);

    % Calculate frequency for each channel
    for i = 1:num_channels
        % Calculate analytic signal using custom Hilbert transform
        analytic_signal = custom_hilbert(voltage_data(:,i));

        % Calculate instantaneous phase
        instantaneous_phase = unwrap(angle(analytic_signal));

        % Calculate instantaneous frequency
        freq_data(:,i) = diff(instantaneous_phase) / (2.0 * pi) / dt;
    end

    % Assign outputs
    frequency.frequency = freq_data;
    frequency.time = time(1:end-1);
end

function analytic = custom_hilbert(x)
    % Compute Hilbert transform without using Signal Processing Toolbox
    N = length(x);
    X = fft(x);

    % Create vector for Hilbert transform in frequency domain
    h = ones(N, 1);
    if mod(N,2) == 0
        % even length
        h(2:N/2) = 2;
        h(N/2+1:end) = 0;
    else
        % odd length
        h(2:(N+1)/2) = 2;
        h((N+1)/2+1:end) = 0;
    end

    % Calculate analytic signal
    analytic = ifft(X .* h);
end
