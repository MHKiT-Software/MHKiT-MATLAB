function figure=plot_current_timeseries(data, principal_direction, ...
                                       varargin)
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
% 
% Returns
% ---------
%   figure: timeseries plot of current-speed velocity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotate coordinate system by supplied principal_direction
principal_directions = data.d - principal_direction;

% Calculate the velocity
velocities = data.s .* cos(pi/180*principal_directions);

% converting epoc time to matlab datetime
time = datetime(data.time, 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS');

% Call on standard xy plotting
figure = plot(time, velocities);
datetick('x',2,'keeplimits');

grid on;
%axis 'tight';
xlabel('Time [date]','FontSize',20);
ylabel('Velocity [\itm/s\rm]','FontSize',20);

if nargin == 3
    title(varargin{1})
end