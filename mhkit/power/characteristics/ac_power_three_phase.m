function power_ac = ac_power_three_phase(voltage, current, power_factor, varargin)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculates three-phase active (real) power [W] from instantaneous
%     voltage and current measurements. Valid for both balanced and
%     unbalanced three-phase systems.
%
% Parameters
% ------------
%   voltage: structure
%       voltage.voltage : Matrix of voltage measurements [V] (3 phases)
%       voltage.time : Time vector
%   current: structure
%       current.current : Matrix of current measurements [A] (3 phases)
%       current.time : Time vector
%   power_factor: double
%       Power factor for the efficiency of the system (0 to 1)
%   Name-Value Parameters:
%       'LineToLine': logical (default: false)
%           Set to true if the given voltage measurements are line-to-line
%       'VoltageImbalanceThreshold': double (default: 2.0)
%           Maximum allowable voltage unbalance factor in percent
%           Default per IEC TS 62600-30:2018, Section 7.1.3
%           Can be overridden based on local system operator's codes
%
% Returns
% ---------
%   power_ac: structure
%       power_ac.power : Vector of calculated active power [W]
%       power_ac.time : Time vector
%
% Standards Reference
% ------------------
% IEC TS 62600-30:2018 Section 7.1.3 Test conditions states:
% "The voltage unbalance factor shall be less than a value, defined by the
% local system operator's codes and standards, measured as 10 min data at
% the MEC unit terminals, or PCC as appropriate. Where no such codes or
% standards exist, a value of 2 % shall be used. The voltage unbalance
% factor may be determined as described in IEC 61800-3:2017, Clause B.5."
%
% Key Equations
% -------------
% 1. Line-to-Line to Phase Voltage Conversion:
%    V_phase = V_line-to-line / √3
%
% 2. Three-Phase Active Power:
%    P = (V1×I1 + V2×I2 + V3×I3) × power_factor
%    where V1,V2,V3 are instantaneous phase voltages
%    and I1,I2,I3 are instantaneous phase currents
%
% 3. Voltage Unbalance Factor:
%    Unbalance = (max_value - min_value) / mean_value × 100%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse input parameters
    p = inputParser;
    addParameter(p, 'LineToLine', false, @islogical);
    addParameter(p, 'VoltageImbalanceThreshold', 2.0, @(x) isnumeric(x) && x > 0);
    parse(p, varargin{:});

    % Extract parameters
    line_to_line = p.Results.LineToLine;
    v_threshold = p.Results.VoltageImbalanceThreshold;

    % Input validation
    if ~isstruct(voltage) || ~isstruct(current)
        error('voltage and current must be structures with voltage/current and time fields');
    end

    % Extract data
    v_data = voltage.voltage;
    i_data = current.current;
    time = voltage.time;

    % Validate data types
    if ~isnumeric(v_data) || ~isnumeric(i_data)
        error('Voltage and current data must be numeric');
    end
    if ~isnumeric(power_factor)
        error('Power factor must be numeric');
    end

    % Validate power factor range
    if power_factor < 0 || power_factor > 1
        error('Power factor must be between 0 and 1');
    end

    % Check for empty data
    if isempty(v_data) || isempty(i_data)
        error('Voltage and current data cannot be empty');
    end

    % Verify three phases
    if size(v_data, 2) ~= 3 || size(i_data, 2) ~= 3
        error('voltage and current must have exactly three phases (columns)');
    end

    % Verify matching dimensions
    if size(v_data, 1) ~= size(i_data, 1)
        error('Voltage and current measurements must have the same number of samples');
    end

    % Verify time vectors match
    if ~isequal(voltage.time, current.time)
        error('Time vectors in voltage and current structures must match');
    end

    % Check for voltage unbalance factor
    v_mean = mean(v_data, 2);
    v_max = max(v_data, [], 2);
    v_min = min(v_data, [], 2);
    v_unbalance = (v_max - v_min) ./ v_mean * 100;

    if any(v_unbalance > v_threshold)
        warning(['Voltage unbalance factor exceeds %.1f%% (max: %.1f%%)\n' ...
                'IEC TS 62600-30:2018 requires voltage unbalance factor ' ...
                'to be below local system operator defined threshold, ' ...
                'or 2%% where no such threshold exists.'], ...
                v_threshold, max(v_unbalance));
    end

    % Convert line-to-line voltage to phase voltage if necessary
    if line_to_line
        v_data = v_data / sqrt(3);
    end

    % Calculate instantaneous power for each phase
    power = sum(v_data .* i_data, 2);

    % Apply power factor
    power = abs(power) * power_factor;

    % Assign outputs
    power_ac = struct();
    power_ac.power = power;
    power_ac.time = time;
end
