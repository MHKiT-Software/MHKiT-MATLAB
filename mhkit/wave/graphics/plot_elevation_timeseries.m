function figure=plot_elevation_timeseries(wave_elevation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots wave elevation timeseries
%     
% Parameters
% ------------
%     wave_elevation: Structure of the following form:
%
%         wave_elevation.elevation: elevation [m]
%
%         wave_elevation.time: time (s);
%         
% Returns
% ---------
%     figure: figure
%         Plot of wave elevation vs. time
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure=plot(wave_elevation.time,wave_elevation.elevation);
xlabel('Time (s)')
ylabel('Wave Elevation (m)') 

