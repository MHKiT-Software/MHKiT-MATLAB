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

py.importlib.import_module('mhkit');

result=py.mhkit.river.performance.rectangular(h,w);

resultc=cell(result);
D_E=resultc{1};
projected_capture_area=resultc{2};

