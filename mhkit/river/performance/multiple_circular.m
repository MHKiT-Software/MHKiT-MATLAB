function [D_E, projected_capture_area] = multiple_circular(diameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the equivalent diameter and projected capture area of a
% multiple circular turbine
%
% Parameters
% ------------
% diameters : vector [m]
%   Device diameters
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
    diameters {mustBeNumeric, mustBeVector}
end

arguments (Output)
    D_E (1,1) {mustBeNumeric}
    projected_capture_area (1,1) {mustBeNumeric}
end

diameters_squared = diameters.^2;
D_E = sqrt(sum(diameters_squared));
projected_capture_area = 0.25 * pi * sum(diameters_squared);

end
