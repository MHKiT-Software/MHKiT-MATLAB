function [D_E,projected_capture_area]=rectangular(h,w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the equivalent diameter and projected capture area of a
%     retangular turbine
%
% Parameters
% ------------
%     h : float
%         Turbine height [m]
%
%     w : float
%         Turbine width [m]
%
% Returns
% ---------
%     D_E : float
%        Equivalent diameter [m]
%
%     projected_capture_area : float
%         Projected capture area [m^2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    h (1,1) {mustBeNumeric}
    w (1,1) {mustBeNumeric}
end
arguments (Output)
    D_E (1,1) {mustBeNumeric}
    projected_capture_area (1,1) {mustBeNumeric}
end

D_E = sqrt(4.0 * h * w / pi);
projected_capture_area = h * w;

end

