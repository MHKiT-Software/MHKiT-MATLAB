function [ebb,flood]=principal_flow_direction(directions,width_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the principal flow directions of current data
%     The weighted average (over the working velocity range of the TEC) 
%     should be considered to be the principal direction of the current, 
%     and should be used for both the ebb and flood cycles to determine 
%     the TEC optimum orientation. 
%     
%     Parameters
%     ----------
%     directions: vector
%       flow directions [degrees]
%         
%     width_dir: int or vector
%       width of direction bins [degrees]
%     Returns
%     -------
%     ebb: float
%         principal ebb direction [degrees]
%     flood: float
%         principal flood direction [degrees]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

directions = py.numpy.array(directions);
dir = py.mhkit.tidal.resource.principal_flow_direction(directions,width_dir);
dirc=cell(dir);
ebb = dirc{1};
flood = dirc{2};
