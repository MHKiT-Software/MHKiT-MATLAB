function Fr=Froude_number(v,h,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculate the Froude Number of the river, channel or duct flow,
%     to check subcritical flow assumption (if Fr <1).
%
% Parameters
% ------------
%     v : float
%         Average Velocity [m/s].
%
%     h : float
%         Mean hydrolic depth float [m].
%
%     g : float (optional)
%         gravitational acceleration [m/s2].
%
% Returns
% ---------
%     Fr : float
%         Froude Number of the river [unitless].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    v (1,1) {mustBeNumeric}
    h (1,1) {mustBeNumeric}
    g (1,1) {mustBeNumeric} = 9.80665
end
arguments (Output)
    Fr (1,1) {mustBeNumeric}
end

Fr = v / sqrt(g * h);

end

