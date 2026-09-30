function [D_E, projected_capture_area] = rectangular(h, w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the equivalent diameter and projected capture area of a
% rectangular turbine
%
% Parameters
% ------------
% h : double [m]
%   Turbine height
% w : double [m]
%   Turbine width
%
% Returns
% ---------
% D_E : double [m]
%   Equivalent diameter
% projected_capture_area : double [m^2]
%   Projected capture area
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
