function Cp = power_coefficient(power, inflow_speed, capture_area, rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the power coefficient of a MEC device
%
% Parameters
% ------------
% power : numeric [W]
%   Power output signal of device after losses
% inflow_speed : numeric [m/s]
%   Velocity of inflow condition
% capture_area : double [m^2]
%   Projected area of rotor normal to inflow
% rho : double [kg/m^3]
%   Density of environment
%
% Returns
% ---------
% Cp : numeric [-]
%   Power coefficient of device
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
