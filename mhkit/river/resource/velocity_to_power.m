function p=velocity_to_power(V,polynomial_coefficients,cut_in,cut_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates power given velocity data and the relationship
%     between velocity and power from an individual turbine
%
% Parameters
% ----------
%     V : Velocity [m/s]
%
%         structure of form:
%
%           V.V: Velocity [m/s]
%
%           V.time: time [datetime or s]
%
%     polynomial_coefficients : vector
%         vector of polynomial coefficients that discribe the relationship between
%         velocity and power at an individual turbine
%
%     cut_in: float
%         Velocity values below cut_in are not used to compute P
%
%     cut_out: float
%         Velocity values above cut_out are not used to compute P
%
% Returns
% -------
%     p : Structure
%
%
%        P.P: Power [W]
%
%        P.time: epoch time [s]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

