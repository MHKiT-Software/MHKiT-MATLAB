function figure=plot_tidal_phase_probability(data, flood, ebb, ...
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

decimals = round(options.bin_size/0.1);
N_bins = int32(round(max(data.s),decimals)/0.1);

[H, bins] = histcounts(data.s, N_bins);
H_ebb = histcounts(data.s(isEbb), bins);
H_flood = histcounts(data.s(~isEbb), bins);

p_ebb = H_ebb./H;
p_flood = H_flood./H;

center = (bins(1:end-1) + bins(2:end))/2;
p_ebb_less = find(p_ebb<p_flood);
p_ebb_greater = find(p_flood<p_ebb);

figure = bar(center([p_ebb_greater]), p_ebb([p_ebb_greater]), 0.9,"blue");

hold on

bar(center([p_ebb_greater]), p_flood([p_ebb_greater]), 0.9,...
    'FaceColor',[0.9290 0.6940 0.1250]);
bar(center([p_ebb_less]), p_flood([p_ebb_less]), 0.9,...
    'FaceColor',[0.9290 0.6940 0.1250]);
bar(center([p_ebb_less]), p_ebb([p_ebb_less]), 0.9, "blue");
xlabel('Velocity [\itm/s\rm]','FontSize',20);
ylabel('Probability','FontSize',20);
legend('Ebb','Flood')
grid on
title(options.title)

hold off
