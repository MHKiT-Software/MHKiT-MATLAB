function [D_E,projected_capture_area]=circular(diameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     Calculates the equivalent diameter and projected capture area of a 
%     circular turbine
%     
% Parameters
% ------------
%     diameter : float
%         Turbine diameter [m]
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

result=py.mhkit.river.device.circular(diameter);

resultc=cell(result);
D_E=resultc{1};
projected_capture_area=resultc{2};