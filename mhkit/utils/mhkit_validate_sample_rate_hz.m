function result = mhkit_validate_sample_rate_hz(time_vector, expected_sample_rate_hz, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Validate that time vector sample rate matches expected rate within tolerance,
% and check whether the sample interval is uniform across the record.
% Handles various time formats: seconds elapsed, POSIX time, MATLAB datenum, etc.
%
% Parameters
% ----------
%   time_vector: array
%       Time vector in any format (seconds elapsed, POSIX time, MATLAB datenum, etc.)
%       Must be monotonically increasing
%
%   expected_sample_rate_hz: double (optional, default = 1 / mean(dt))
%       Expected sample rate in Hz. When omitted, the nominal rate of the
%       record is used, so only the uniformity check is meaningful.
%
%   options.tolerance_percent: double (optional, default = 1.0)
%       Allowed tolerance as percentage from expected rate and from the
%       mean sample interval (0-100)
%
% Returns
% -------
%   result: structure
%       result.pass : logical - true if within tolerance
%       result.min_sample_rate_hz : double - minimum observed sample rate
%       result.median_sample_rate_hz : double - median observed sample rate (used for pass/fail)
%       result.mean_sample_rate_hz : double - mean observed sample rate
%       result.max_sample_rate_hz : double - maximum observed sample rate
%       result.std_sample_rate_hz : double - standard deviation of the observed sample rate
%       result.deviation_percent : double - deviation of the median sample rate from the expected rate (used for pass/fail)
%       result.max_deviation_percent : double - largest deviation of any single sample interval (diagnostic only)
%       result.is_uniform : logical - true if every sample interval is within tolerance of the mean interval
%       result.max_interval_deviation_percent : double - largest deviation of any sample interval from the mean interval (used for is_uniform)
%       result.tolerance_percent : double - tolerance percentage used
%       result.time_format : string - detected time format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        time_vector double {mustBeVector}
        expected_sample_rate_hz double {mustBePositive} = []
        options.tolerance_percent double {mustBeInRange(options.tolerance_percent, 0, 100)} = 1.0
    end
    
    % Validate minimum length
    if length(time_vector) < 2
        error('MHKiT:mhkit_validate_sample_rate_hz: time_vector must have at least 2 elements');
    end
    
    % Ensure column vector for consistency
    time_vector = time_vector(:);
    
    % Validate monotonically increasing
    if any(diff(time_vector) <= 0)
        error('MHKiT:mhkit_validate_sample_rate_hz: time_vector must be monotonically increasing');
    end
    
    % Detect and handle time format
    [time_seconds, time_format] = detect_and_convert_time_format(time_vector);
    
    % Calculate time differences (dt) in seconds
    time_diffs_s = diff(time_seconds);
    
    % Calculate sample rates from time differences
    % sample_rate = 1 / dt
    sample_rates_hz = 1 ./ time_diffs_s;
    
    % Calculate statistics
    min_sample_rate_hz = min(sample_rates_hz);
    mean_sample_rate_hz = mean(sample_rates_hz);
    median_sample_rate_hz = median(sample_rates_hz);
    max_sample_rate_hz = max(sample_rates_hz);
    std_sample_rate_hz = std(sample_rates_hz);
    
    % Uniformity: every interval within tolerance of the mean interval
    mean_time_diff_s = mean(time_diffs_s);
    max_interval_deviation_percent = max(abs(time_diffs_s - mean_time_diff_s)) / mean_time_diff_s * 100;
    is_uniform = max_interval_deviation_percent <= options.tolerance_percent;
    
    if isempty(expected_sample_rate_hz)
        expected_sample_rate_hz = 1 / mean_time_diff_s;
    end
    
    % Calculate deviations from expected rate
    min_rate_deviation_percent = abs(min_sample_rate_hz - expected_sample_rate_hz) / expected_sample_rate_hz * 100;
    max_rate_deviation_percent = abs(max_sample_rate_hz - expected_sample_rate_hz) / expected_sample_rate_hz * 100;
    mean_rate_deviation_percent = abs(mean_sample_rate_hz - expected_sample_rate_hz) / expected_sample_rate_hz * 100;
    median_rate_deviation_percent = abs(median_sample_rate_hz - expected_sample_rate_hz) / expected_sample_rate_hz * 100;
    
    % Largest deviation of any single sample interval (reported as a diagnostic)
    max_deviation_percent = max([min_rate_deviation_percent, max_rate_deviation_percent, mean_rate_deviation_percent]);
    
    % Pass/fail is judged on the median sample rate, which is tolerant of
    % outliers. This catches a wrong expected_sample_rate_hz argument without
    % failing on the occasional irregular timestamp interval present in real
    % data acquisition systems.
    pass = median_rate_deviation_percent <= options.tolerance_percent;
    
    % Create result structure
    result = struct();
    result.pass = pass;
    result.min_sample_rate_hz = min_sample_rate_hz;
    result.median_sample_rate_hz = median_sample_rate_hz;
    result.mean_sample_rate_hz = mean_sample_rate_hz;
    result.max_sample_rate_hz = max_sample_rate_hz;
    result.std_sample_rate_hz = std_sample_rate_hz;
    result.deviation_percent = median_rate_deviation_percent;
    result.max_deviation_percent = max_deviation_percent;
    result.is_uniform = is_uniform;
    result.max_interval_deviation_percent = max_interval_deviation_percent;
    result.tolerance_percent = options.tolerance_percent;
    result.time_format = time_format;

end

function [time_seconds, time_format] = detect_and_convert_time_format(time_vector)
    % Detect time format and convert to seconds elapsed
    
    % Get time range and typical values
    time_range = time_vector(end) - time_vector(1);
    first_value = time_vector(1);
    
    % Detection logic based on typical ranges and values
    if first_value > 1e9
        % POSIX timestamp (seconds since 1970-01-01)
        % Typical range: 1e9 to 2e9 for years 2001-2033
        time_format = "POSIX timestamp";
        time_seconds = time_vector - time_vector(1); % Convert to elapsed seconds
        
    elseif first_value > 1e5 && first_value < 1e9
        % Likely MATLAB datenum (days since year 0000)
        % Typical range: ~700,000 to ~800,000 for modern dates
        time_format = "MATLAB datenum";
        time_seconds = (time_vector - time_vector(1)) * 24 * 3600; % Convert days to seconds
        
    elseif first_value >= 0 && time_range < 1e6
        % Likely seconds elapsed (starting from 0 or small value)
        time_format = "seconds elapsed";
        time_seconds = time_vector - time_vector(1); % Normalize to start at 0
        
    else
        % Unknown format - assume it's already in seconds and warn
        time_format = "unknown (assuming seconds)";
        time_seconds = time_vector - time_vector(1);
        warning('MHKiT:mhkit_validate_sample_rate_hz: Unknown time format, assuming seconds elapsed');
    end
    
    % Final validation: ensure we have positive time differences
    if any(diff(time_seconds) <= 0)
        error('MHKiT:mhkit_validate_sample_rate_hz: Time conversion resulted in non-increasing values');
    end
    
end
