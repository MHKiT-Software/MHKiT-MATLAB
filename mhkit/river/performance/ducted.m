function [D_E,projected_capture_area]=ducted(diameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the equivalent diameter and projected capture area of a
%     ducted turbine
%
% Parameters
% ------------
%     diameter : float
%         ducted diameter [m]
%
% Returns
% ---------
%     D_E : float
%        Equivalent diameter [m]
%
%     projected_capture_area : float
%         Projected capture area [m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');

result=py.mhkit.river.performance.ducted(diameter);

resultc=cell(result);
D_E=resultc{1};
projected_capture_area=resultc{2};

