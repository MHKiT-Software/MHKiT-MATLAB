function [M_flap,M_edge] = blade_moments(blade_coefficients,flap_offset,...
    flap_raw, edge_offset, edge_raw)

%%%%%%%%%%%%%%%%%%%%
%     Transfer function for deriving blade flap and edge moments using blade matrix.
%
% Parameters
% ------------
%     blade_coefficients: vector
%         Derived blade calibration coefficients listed in order of D1, D2, D3, D4
%
%     flap_offset : double or int
%         Derived offset of raw flap signal obtained during calibration process
%
%     flap_raw : vector
%         Raw strain signal of blade in the flapwise direction
%
%     edge_offset : double or int
%         Derived offset of raw edge signal obtained during calibration process
%
%     edge_raw : vector
%         Raw strain signal of blade in the edgewise direction
%
% Returns
% ---------
%     M_flap: vector
%         Blade flapwise moment in SI units
%
%     M_edge: vector
%         Blade edgewise moment in SI units
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

data = py.mhkit.loads.general.blade_moments(py.numpy.array(blade_coefficients),...
    flap_offset,py.numpy.array(flap_raw),edge_offset,py.numpy.array(edge_raw));

data_cell = cell(data);

M_flap = double(py.array.array('d',py.numpy.nditer(data_cell{1},pyargs("flags",{"refs_ok"}))));

M_edge = double(py.array.array('d',py.numpy.nditer(data_cell{2},pyargs("flags",{"refs_ok"}))));

