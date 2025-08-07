function frequency = instantaneous_frequency(voltage)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates instantaneous frequency of measured voltage
%     
% Parameters
% ------------
%   voltage: structure
%       voltage.voltage : Measured voltage data (V) with each timeseries in its own column
%       voltage.time : Time vector corresponding to voltage measurements
%
% Returns
% ---------
%   frequency: structure
%       frequency.frequency : Instantaneous frequency of the measured voltage (Hz)
%       frequency.time : Time vector (one element shorter than input due to differentiation)
%
% Key Equations
% -------------
% 1. Analytic signal using Hilbert transform:
%    z(t) = voltage(t) + i * H[voltage(t)]
%
% 2. Instantaneous phase:
%    φ(t) = unwrap(angle(z(t)))
%
% 3. Instantaneous frequency:
%    f(t) = (1/(2π)) * dφ/dt
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStruct = @(x) isstruct(x);
    
    % Add required parameters
    addRequired(p, 'voltage', validStruct);
    
    % Parse inputs
    parse(p, voltage);
    
    % Extract validated inputs
    voltage = p.Results.voltage;
    
    % Validate input structure has required fields
    if ~isfield(voltage, 'voltage')
        error('MHKiT:instantaneous_frequency: voltage structure must contain voltage field');
    end
    if ~isfield(voltage, 'time')
        error('MHKiT:instantaneous_frequency: voltage structure must contain time field');
    end
    
    % Extract data from structure
    voltage_data = voltage.voltage;
    time_vector = voltage.time;
    
    % Validate dimensions
    if size(voltage_data, 1) ~= length(time_vector)
        error('MHKiT:instantaneous_frequency: voltage data rows must match time vector length');
    end
    
    % Get data dimensions
    [num_samples, num_columns] = size(voltage_data);
    
    % Validate minimum data length for meaningful frequency calculation
    if num_samples < 4
        error('MHKiT:instantaneous_frequency: voltage data must have at least 4 samples for frequency calculation');
    end
    
    % Calculate time differences for frequency calculation
    time_diff = diff(time_vector(:));  % Ensure column vector
    
    % Check for uniform time spacing (within tolerance)
    dt_mean = mean(time_diff);
    % Warn the user if the sample rate is not consistent
    sample_rate_hz_tolerance_percent = 0.01;  % 1% tolerance
    dt_tolerance = 0.01 * dt_mean;  % 1% tolerance
    if any(abs(time_diff - dt_mean) > dt_tolerance)
        warning('MHKiT:instantaneous_frequency:SampleRateVariation', ...
                ['The sample rate of this signal has differences greater than ', ...
                 num2str(sample_rate_hz_tolerance_percent * 100), '%% from the mean sample rate. ', ...
                 'This has a high likelihood of causing inaccurate instantaneous frequency calculation results']);
        warning('MHKiT:instantaneous_frequency:MeanSampleRate', ...
                ['Mean sample rate hz: ', num2str(1/dt_mean)]);
        warning('MHKiT:instantaneous_frequency:MaxSampleRate', ...
                ['Max sample rate hz: ', num2str(1/min(time_diff))]);
        warning('MHKiT:instantaneous_frequency:MinSampleRate', ...
                ['Min sample rate hz: ', num2str(1/max(time_diff))]);
        warning('MHKiT:instantaneous_frequency:StdSampleRate', ...
                ['Standard Deviation sample rate hz: ', num2str(std(1./time_diff))]);
    end
    
    % Initialize output frequency matrix
    frequency_data = zeros(num_samples - 1, num_columns);
    
    % Process each column of voltage data
    for col_idx = 1:num_columns
        current_voltage = voltage_data(:, col_idx);
        
        % Check if Signal Processing Toolbox is available for hilbert function
        if ~exist('hilbert', 'file')
            warning('MHKiT:instantaneous_frequency: MATLAB hilbert function not available, using scipy hilbert function');
            analytic_signal = mhkit_scipy_hilbert(current_voltage);
        else
            % Apply Hilbert transform to get analytic signal
            analytic_signal = hilbert(current_voltage);
        end
        
        % Calculate instantaneous phase
        instantaneous_phase = angle(analytic_signal);
        
        % Unwrap phase to remove 2π discontinuities
        unwrapped_phase = unwrap(instantaneous_phase);
        
        % Calculate instantaneous frequency
        phase_diff = diff(unwrapped_phase);
        % Ensure proper element-wise division with matching dimensions
        if size(phase_diff, 1) == 1
            phase_diff = phase_diff(:);  % Convert to column vector if row
        end
        instantaneous_frequency = phase_diff ./ (2.0 * pi * time_diff);
        
        % Store result
        frequency_data(:, col_idx) = instantaneous_frequency;
    end
    
    % Create output structure
    frequency = struct();
    frequency.frequency = frequency_data;
    frequency.time = time_vector(2:end);  % Time vector is one element shorter due to differentiation

end
