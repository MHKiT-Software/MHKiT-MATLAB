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

py.importlib.import_module('mhkit');
diameters=py.list(diameters);
result=py.mhkit.river.performance.multiple_circular(diameters);

resultc=cell(result);
D_E=resultc{1};
projected_capture_area=resultc{2};

