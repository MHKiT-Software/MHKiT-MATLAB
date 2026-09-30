function P = ac_power_three_phase(voltage, current, power_factor, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate magnitude of active AC power from line to neutral voltage and current
%
% Parameters
% ------------
%   voltage: structure or matrix
%       voltage.voltage : Three-phase voltage measurements [V] (n_time x 3 matrix)
%                        Each row represents one time step, columns are phases A, B, C
%       voltage.time : Time vector (n_time x 1) (if time series data)
%   current: structure or matrix
%       current.current : Three-phase current measurements [A] (n_time x 3 matrix)
%                        Each row represents one time step, columns are phases A, B, C
%       current.time : Time vector (n_time x 1) (if time series data)
%   power_factor: numeric scalar
%       Power factor for the efficiency of the system [dimensionless]
%   'LineToLine': name-value pair (optional)
%       Set to true if voltage measurements are line-to-line (default: false)
%
% Returns
% ---------
%   P: structure
%       P.power : Magnitude of active AC power [W]
%       P.time : Time vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStructOrNumeric = @(x) isstruct(x) || isnumeric(x);
    validNumeric = @(x) isnumeric(x) && isscalar(x);
    validLogical = @(x) islogical(x) && isscalar(x);
    
    % Add required parameters
    addRequired(p, 'voltage', validStructOrNumeric);
    addRequired(p, 'current', validStructOrNumeric);
    addRequired(p, 'power_factor', validNumeric);
    
    % Add optional name-value pairs
    addParameter(p, 'LineToLine', false, validLogical);
    
    % Parse inputs
    parse(p, voltage, current, power_factor, varargin{:});
    
    % Extract validated inputs
    voltage = p.Results.voltage;
    current = p.Results.current;
    power_factor = p.Results.power_factor;
    line_to_line = p.Results.LineToLine;

    % Validate power factor is between 0 and 1
    if power_factor < 0 || power_factor > 1
        error('MHKiT:ac_power_three_phase: power_factor must be between 0 and 1 (inclusive). Received: %.3f', power_factor);
    end

    % Extract data and time vectors
    if isstruct(voltage)
        % Validate input structures have required fields
        if ~isfield(voltage, 'voltage')
            error('MHKiT:ac_power_three_phase: voltage structure must contain voltage field');
        end
        if ~isfield(voltage, 'time')
            error('MHKiT:ac_power_three_phase: voltage structure must contain time field');
        end
        voltage_data = voltage.voltage;
        voltage_time = voltage.time;
    else
        voltage_data = voltage;
        voltage_time = (1:size(voltage, 1))';  % Default time vector
    end
    
    if isstruct(current)
        % Validate input structures have required fields
        if ~isfield(current, 'current')
            error('MHKiT:ac_power_three_phase: current structure must contain current field');
        end
        if ~isfield(current, 'time')
            error('MHKiT:ac_power_three_phase: current structure must contain time field');
        end
        current_data = current.current;
        current_time = current.time;
    else
        current_data = current;
        current_time = (1:size(current, 1))';  % Default time vector
    end

    % Validate voltage has three columns
    if size(voltage_data, 2) ~= 3
        error('MHKiT:ac_power_three_phase: voltage must have three columns for three-phase measurements');
    end
    
    % Validate current has three columns
    if size(current_data, 2) ~= 3
        error('MHKiT:ac_power_three_phase: current must have three columns for three-phase measurements');
    end
    
    % Validate dimensions match
    if ~isequal(size(voltage_data), size(current_data))
        error('MHKiT:ac_power_three_phase: voltage and current must have the same dimensions');
    end

    % Validate time vectors match
    if ~isequal(voltage_time, current_time)
        error('MHKiT:ac_power_three_phase: voltage and current time vectors must match');
    end

    % After validation, use the time vector from voltage
    time_vector = voltage_time;

    % Instantaneous power summed across the three phases, as in MHKiT-Python
    active_power = abs(sum(voltage_data .* current_data, 2)) * power_factor;

    if line_to_line
        active_power = active_power * sqrt(3);
    end

    % Create output structure
    P = struct();
    P.power = active_power;
    P.time = time_vector;

end
