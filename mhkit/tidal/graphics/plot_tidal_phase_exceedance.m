function figure=plot_tidal_phase_exceedance(data, flood, ebb, ...
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
%    flood: float
%        principal flood direction [degrees]
%
%    ebb: float
%        principal ebb direction [degrees]
% 
%    bin_size: numeric (optional)
%       Speed bin size. Default = 0.1 m/s
%       to call: plot_tidal_phase_probability(data, flood, ebb,"bin_size",bin_size)
%
%    title: string (optional)
%       title for the plot 
%       to call: plot_tidal_phase_probability(data, flood, ebb,"title",title)
%
% Returns
% ---------
%   figure: stacked bar graph of the probability of exceedance in 
%           flood and ebb directions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data
    flood
	ebb
    options.bin_size = 0.1; % m/s
    options.title = "";
end

%check to see if the first input argument is a structure
if any(~isstruct(data))
    ME = MException('MATLAB:plot_tidal_phase_probability','data must be a structure');
    throw(ME);
end

%check to see if the second input argument is a number
if any([~isnumeric(flood), length(flood) ~= 1])
    ME = MException('MATLAB:plot_tidal_phase_probability','flood must be a number');
    throw(ME);
end

%check to see if the third input argument is a number
if any([~isnumeric(ebb), length(ebb) ~= 1])
    ME = MException('MATLAB:plot_tidal_phase_probability','ebb must be a number');
    throw(ME);
end

max_angle = max(ebb, flood);
min_angle = min(ebb, flood);
    
lower_split = rem((min_angle + (360 - max_angle + min_angle)/2 ) , 360);
upper_split = lower_split + 180;
    
if (lower_split <= ebb) && (ebb < upper_split)
	isEbb = ((data.d < upper_split) & (data.d >= lower_split));
else
	isEbb = ~((data.d < upper_split) & (data.d >= lower_split));
end

s_ebb = struct('time', {data.time(isEbb)},...
    'Discharge', {data.Discharge(isEbb)});
s_flood = struct('time', {data.time(~isEbb)},...
    'Discharge', {data.Discharge(~isEbb)});

F = exceedance_probability(data).F;
F_ebb = exceedance_probability(s_ebb).F;
F_flood = exceedance_probability(s_flood).F;

decimals = round(options.bin_size/0.1);
s_new = [round(min(data.s),decimals):options.bin_size:...
    round(max(data.s),decimals)+options.bin_size];
s_new = s_new(1:end-1); % this accounts for the difference between matlab 
                        % and numpy.arange which uses a half open interval

py.importlib.import_module('scipy');
kwa = pyargs('bounds_error',false);

f_ebb = py.scipy.interpolate.interpolate.interp1d(s_ebb.Discharge,...
    F_ebb, kwa);
f_flood = py.scipy.interpolate.interpolate.interp1d(s_flood.Discharge,...
    F_flood, kwa);

F_ebb = double(f_ebb(s_new));
F_flood = double(f_flood(s_new));

py.importlib.import_module('numpy');
F_max_total = py.numpy.nanmax(F_ebb) + py.numpy.nanmax(F_flood);

s_plot = repmat(s_new,2,1);
figure = area(s_plot.', [F_ebb/F_max_total*100; F_flood/F_max_total*100].');

hold on

xlabel('Velocity [\itm/s\rm]','FontSize',20);
ylabel('Probability of Exceedance','FontSize',18);
legend('Ebb','Flood')
grid on
title(options.title)

hold off
