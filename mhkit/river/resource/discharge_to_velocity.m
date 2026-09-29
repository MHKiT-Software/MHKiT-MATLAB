function V=discharge_to_velocity(Q,polynomial_coefficients)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates velocity given discharge data and the relationship between
%     discharge and velocity at an individual turbine
%
% Parameters
% ------------
%     Q : Discharge data [m3/s]
%
%         structure of form:
%
%            Q.Discharge
%
%            Q.time
%
%     polynomial_coefficients : vector
%         Vector of polynomial coefficients (highest degree first) that
%         describe the relationship between discharge and velocity at an
%         individual turbine
%
% Returns
% ------------
%     V: Structure
%
%
%         V.V: Velocity [m/s]
%
%         V.time: time [datetime or s]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

