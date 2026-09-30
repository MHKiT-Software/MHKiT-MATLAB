function [D_E, projected_capture_area] = ducted(diameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the equivalent diameter and projected capture area of a
% ducted turbine
%
% Parameters
% ------------
% diameter : double [m]
%   Duct diameter
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
    diameter (1,1) {mustBeNumeric}
end

arguments (Output)
    D_E (1,1) {mustBeNumeric}
    projected_capture_area (1,1) {mustBeNumeric}
end

D_E = diameter;
projected_capture_area = (1/4) * pi * (D_E.^2);

end
