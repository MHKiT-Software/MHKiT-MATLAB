function V = discharge_to_velocity(Q, polynomial_coefficients)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates velocity given discharge data and the relationship between
% discharge and velocity at an individual turbine
%
% Parameters
% ------------
% Q : struct
%   Discharge data
%     Q.Discharge : vector or matrix [m^3/s]
%       Discharge
%     Q.time : vector [datetime or s]
%       Time
% polynomial_coefficients : vector
%   Polynomial coefficients (highest degree first, e.g. poly.coef from
%   polynomial_fit) that describe the relationship between discharge and
%   velocity at an individual turbine
%
% Returns
% ---------
% V : struct
%   Velocity data
%     V.V : vector or matrix [m/s]
%       Velocity, one value per discharge value
%     V.time : vector [s]
%       Time, with datetime converted to epoch seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Q struct
    polynomial_coefficients {mustBeNumeric, mustBeVector}
end

arguments (Output)
    V struct
end

time = Q.time;
if any(isdatetime(time))
    time = posixtime(time);
end

V.V = polyval(polynomial_coefficients, Q.Discharge);
V.time = time;

end
