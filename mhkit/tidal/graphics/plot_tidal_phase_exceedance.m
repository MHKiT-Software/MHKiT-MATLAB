function figure=plot_tidal_phase_exceedance(data, flood, ebb, ...
                                                          options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns a stacked area plot of the exceedance probability for the 
%     flood and ebb tidal phases. 
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
%    savepath: string (optional)
%       path and filename to save figure.
%       to call: plot_tidal_phase_probability(data, flood, ebb,"savepath",savepath)
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
    options.savepath = "";
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

s_ebb = data.s(isEbb);
s_flood = data.s(~isEbb);

F_ebb = exceedance_probability(s_ebb);
F_flood = exceedance_probability(s_flood);

decimals = round(options.bin_size/0.1);
s_new = [round(min(data.s),decimals):options.bin_size:...
    round(max(data.s),decimals)+options.bin_size];
s_new = s_new(1:end-1); % this accounts for the difference between matlab 
                        % and numpy.arange which uses a half open interval

[~,i_ebb] = unique(s_ebb);
[~,i_flood] = unique(s_flood);

f_ebb = py.scipy.interpolate.interpolate.interp1d(s_ebb.Discharge,...
    F_ebb, kwa);
f_flood = py.scipy.interpolate.interpolate.interp1d(s_flood.Discharge,...
    F_flood, kwa);

F_max_total = max(F_ebb(~isnan(F_ebb))) + max(F_flood(~isnan(F_flood)));

s_plot = repmat(s_new,2,1);
figure = area(s_plot.', [F_ebb/F_max_total*100; F_flood/F_max_total*100].');

hold on

xlabel('Velocity [\itm/s\rm]','FontSize',20);
ylabel('Probability of Exceedance','FontSize',18);
legend('Ebb','Flood')
grid on
title(options.title)

len = strlength(options.savepath);
if len > 1
    exportgraphics(gca, options.savepath);
end 

hold off