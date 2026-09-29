function TSR=tip_speed_ratio(rotor_speed, rotor_diameter, inflow_speed)

%%%%%%%%%%%%%%%%%%%%
%     Function used to calculate the tip speed ratio (TSR) of a MEC device with rotor
%
%
% Parameters
% ------------
%     rotor_speed : vector
%         Rotor Speed [rps]
%
%     rotor_diameter : double or int
%         diameter -f rotor [m]
%
%     inflow_speed : vector
%         Velocity of inflow condition [m/s]
%
% Returns
% ---------
%     TSR: vector
%         Calculated tip speed ratio (TSR)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

