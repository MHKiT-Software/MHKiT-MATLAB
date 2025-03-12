function power_dc = dc_power(voltage, current)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculates DC power from voltage and current measurements across
%     multiple channels. Supports both single and multi-channel
%     configurations.
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

    % Input validation
    if ~isstruct(voltage) || ~isstruct(current)
        error('voltage and current must be structures with voltage/current and time fields');
    end

    % Validate required fields exist
    required_fields_v = {'voltage', 'time'};
    required_fields_c = {'current', 'time'};

    if ~all(isfield(voltage, required_fields_v))
        error('voltage structure must contain fields: voltage and time');
    end

    if ~all(isfield(current, required_fields_c))
        error('current structure must contain fields: current and time');
    end

    % Extract data
    v_data = voltage.voltage;
    i_data = current.current;
    v_time = voltage.time;
    c_time = current.time;

    % Verify time vectors match
    if ~isequal(v_time, c_time)
        error('Time vectors in voltage and current structures must match');
    end

    % Verify dimensions match
    if ~isequal(size(v_data), size(i_data))
        error('voltage and current must have the same dimensions');
    end

    % Calculate power for each channel
    power = v_data .* i_data;

    % Calculate gross power (sum across channels - dimension 2)
    gross_power = sum(power, 2);

    % Assign outputs
    power_dc = struct();
    power_dc.power = power;
    power_dc.gross = gross_power;
    power_dc.time = v_time;
end
