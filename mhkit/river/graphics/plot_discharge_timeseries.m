function figure=plot_discharge_timeseries(Q,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plots discharge vs time
%
% Parameters
% ------------
%     Q: structure
%
%      Q.Discharge: Discharge [m/s]
%
%      Q.time: epoch time [s]
%
%     title: string (optional)
%       title for the plot
%       to call: plot_discharge_timeseries(Q,"title",title)
%
%     savepath: string (optional)
%         path and filename to save figure.
%         to call: plot_discharge_timeseries(Q,"savepath",savepath)
%
% Returns
% ---------
%     figure: Plot of discharge vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Q
    options.title = "";
    options.savepath = "";
end

figure=plot(datetime(Q.time, 'convertfrom','posixtime'),Q.Discharge);
datetick('x',1,'keeplimits');
grid on
xlabel('Date','FontSize',20)
ylabel('Discharge [m^{3}/s]','FontSize',20)

title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end

hold off

