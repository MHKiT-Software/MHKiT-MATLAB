function figure=plot_elevation_timeseries(wave_elevation, options)

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
%     title: string (optional)
%         title for the plot
%         to call: plot_elevation_timeseries(wave_elevation,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_elevation_timeseries(wave_elevation,"savepath",savepath)
%         
% Returns
% ---------
%     figure: figure
%         Plot of wave elevation vs. time
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    wave_elevation 
    options.title = "";
    options.savepath = "";
end

figure=plot(wave_elevation.time,wave_elevation.elevation);
xlabel('Time (s)')
ylabel('Wave Elevation (m)') 

title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

