function Cp=power_coefficient(power, inflow_speed, capture_area, rho)

%%%%%%%%%%%%%%%%%%%%
%     Function that calculates the power coefficient of MEC device
%
%
% Parameters
% ------------
%     power : vector
%         Power output signal of device after losses [W]
%
%     inflow_speed : vector
%         Velocity of inflow condition [m/s]
%
%     capture_area : double or int
%         Projected area of rotor normal to inflow [m^2]
%
%     rho : double or int
%         Density of environment [kg/m^3]
%
% Returns
% ---------
%     Cp: vector
%         Power coefficient of device [-]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    power {mustBeNumeric}
    inflow_speed {mustBeNumeric}
    capture_area (1,1) {mustBeNumeric}
    rho (1,1) {mustBeNumeric}
end
arguments (Output)
    Cp {mustBeNumeric}
end

% Predicted power from inflow
power_in = 0.5 .* rho .* capture_area .* inflow_speed.^3;

Cp = power ./ power_in;

end

