function TSR = tip_speed_ratio(rotor_speed, rotor_diameter, inflow_speed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the tip speed ratio (TSR) of a MEC device with rotor
%
% Parameters
% ------------
% rotor_speed : numeric [rev/s]
%   Rotor speed
% rotor_diameter : double [m]
%   Diameter of rotor
% inflow_speed : numeric [m/s]
%   Velocity of inflow condition
%
% Returns
% ---------
% TSR : numeric [-]
%   Calculated tip speed ratio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    rotor_speed {mustBeNumeric}
    rotor_diameter (1,1) {mustBeNumeric}
    inflow_speed {mustBeNumeric}
end

arguments (Output)
    TSR {mustBeNumeric}
end

rotor_velocity = rotor_speed .* pi .* rotor_diameter;

TSR = rotor_velocity ./ inflow_speed;

end
