function [D_E,projected_capture_area]=multiple_circular(diameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the equivalent diameter and projected capture area of a
%     multiple circular turbine
%
% Parameters
% ------------
%     diameters: array or vector
%         vector of device diameters [m]
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

