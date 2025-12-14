function power_dc = dc_power(voltage, current)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates DC power from voltage and current measurements
%
% Parameters
% ------------
%   voltage: structure
%       voltage.voltage : Matrix of voltage measurements [V]
%                         Each row is a time point, each column is a channel
%       voltage.time : Time vector
%   current: structure
%       current.current : Matrix of current measurements [A]
%                         Each row is a time point, each column is a channel
%       current.time : Time vector
%
% Returns
% ---------
%   power_dc: structure
%       power_dc.power : Matrix of calculated power for each channel [W]
%                        Same dimensions as input voltage/current
%       power_dc.gross : Vector of gross power (sum of all channels) [W]
%       power_dc.time : Time vector
%
% Key Equations
% -------------
% 1. Channel Power:
%    P_channel = V_channel × I_channel
%
% 2. Gross Power:
%    P_gross = sum(P_channel) across all channels
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    input_parser = inputParser;
    
    % Define validation functions
    valid_struct = @(x) isstruct(x);
    
    % Add required parameters
    addRequired(input_parser, 'voltage', valid_struct);
    addRequired(input_parser, 'current', valid_struct);
    
    % Parse inputs
    parse(input_parser, voltage, current);
    
    % Extract validated inputs
    voltage = input_parser.Results.voltage;
    current = input_parser.Results.current;
    
    % Validate input structures have required fields
    if ~isfield(voltage, 'voltage')
        error('MHKiT:dc_power: voltage structure must contain voltage field');
    end
    if ~isfield(voltage, 'time')
        error('MHKiT:dc_power: voltage structure must contain time field');
    end
    if ~isfield(current, 'current')
        error('MHKiT:dc_power: current structure must contain current field');
    end
    if ~isfield(current, 'time')
        error('MHKiT:dc_power: current structure must contain time field');
    end
    
    % Extract data matrices
    voltage_data = voltage.voltage;
    current_data = current.current;
    voltage_time = voltage.time;
    current_time = current.time;
    
    % Validate dimensions match
    if ~isequal(size(voltage_data), size(current_data))
        error('MHKiT:dc_power: voltage and current must have the same dimensions');
    end
    
    % Validate time vectors match
    if ~isequal(voltage_time, current_time)
        error('MHKiT:dc_power: Time vectors must match between voltage and current structures');
    end
    
    % Calculate power for each channel (element-wise multiplication)
    power_matrix = voltage_data .* current_data;
    
    % Calculate gross power (sum across channels - dimension 2)
    % Handle NaN values by skipping them in the sum
    gross_power_vector = sum(power_matrix, 2, 'omitnan');
    
    % Create output structure
    power_dc = struct();
    power_dc.power = power_matrix;
    power_dc.gross = gross_power_vector;
    power_dc.time = voltage_time;

end
