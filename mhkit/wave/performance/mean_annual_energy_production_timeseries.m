function maep=mean_annual_energy_production_timeseries(L,J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Calculates mean annual energy production (MAEP) from timeseries
%     
% Parameters
% ------------
%     L: numpy array or vector
%         Capture length
%
%     J: numpy array or vector
%         Wave energy flux
%         
% Returns
% ---------
%     maep: float
%         Mean annual energy production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');

J=py.numpy.array(J);
L=py.numpy.array(L);

maep=double(py.mhkit.wave.performance.mean_annual_energy_production_timeseries(L,J));

