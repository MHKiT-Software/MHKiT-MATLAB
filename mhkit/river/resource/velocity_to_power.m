function p = velocity_to_power(V, polynomial_coefficients, cut_in, cut_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates power given velocity data and the relationship between
% velocity and power from an individual turbine
%
% Parameters
% ------------
% V : struct
%   Velocity data
%     V.V : vector or matrix [m/s]
%       Velocity
%     V.time : vector [datetime or s]
%       Time
% polynomial_coefficients : vector
%   Polynomial coefficients (highest degree first, e.g. poly.coef from
%   polynomial_fit) that describe the relationship between velocity and
%   power at an individual turbine
% cut_in : double [m/s]
%   Velocity values below cut_in produce 0 power
% cut_out : double [m/s]
%   Velocity values above cut_out produce 0 power
%
% Returns
% ---------
% p : struct
%   Power data
%     p.P : vector or matrix [W]
%       Power, one value per velocity value
%     p.time : vector [s]
%       Time, with datetime converted to epoch seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    V struct
    polynomial_coefficients {mustBeNumeric, mustBeVector}
    cut_in (1,1) {mustBeNumeric}
    cut_out (1,1) {mustBeNumeric}
end

arguments (Output)
    p struct
end

time = V.time;
if any(isdatetime(time))
    time = posixtime(time);
end

velocity = V.V;
power = polyval(polynomial_coefficients, velocity);

% Turbine produces 0 power outside of the cut-in/cut-out bounds
power(velocity < cut_in) = 0.0;
power(velocity > cut_out) = 0.0;

p.P = power;
p.time = time;

end
