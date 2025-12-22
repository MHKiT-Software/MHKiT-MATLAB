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


% Flatten arrays to 1D to handle both row and column vectors
J=py.numpy.array(J).flatten();
L=py.numpy.array(L).flatten();

maep=double(py.mhkit.wave.performance.mean_annual_energy_production_timeseries(L,J));

