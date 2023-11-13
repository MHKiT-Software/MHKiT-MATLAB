function figure=plot_current_timeseries(data, principal_direction, ...
                                       options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns a plot of velocity from an array of direction and speed
%     data in the direction of the supplied principal_direction.
%
% Parameters
% ------------
%    data: structure
%
%      data.time: vector
%       days from January 0, 0000 in the proleptic ISO calendar
%
%      data.d: vector
%       time-series of directions [degrees]
%
%      data.s: vector
%       time-series of speeds [cm/s]
%
%    principal_direction: numeric
%        Direction to compute the velocity in [degrees]
%
%    title: string (optional)
%       title for the plot
%       to call: plot_current_timeseries(data,principal_direction,"title",title)
%
%    savepath: string (optional)
%       path and filename to save figure.
%       to call: plot_current_timeseries(data,principal_direction,"savepath",savepath)
%
% Returns
% ---------
%   figure: timeseries plot of current-speed velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data
    principal_direction
    options.title = "";
    options.savepath = "";
end

% Rotate coordinate system by supplied principal_direction
principal_directions = data.d - principal_direction;

% Calculate the velocity
velocities = data.s .* cos(pi/180*principal_directions);

% converting epoc time to matlab datetime
time = datetime(data.time, 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS');

% Call on standard xy plotting
figure = plot(time, velocities);
hold on
datetick('x',2,'keeplimits');

grid on;
%axis 'tight';
xlabel('Time [date]','FontSize',20);
ylabel('Velocity [\itm/s\rm]','FontSize',20);

title(options.title)

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end

hold off

