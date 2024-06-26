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

directions = data.d; % degrees
timestamps = data.time;
velocities = data.s; % cm per second

    isEbb = flood_or_ebb(directions, flood, ebb);

    ebb_velocities = velocities(isEbb);
    ebb_timestamps = timestamps(isEbb);

    flood_velocities = velocities(~isEbb);
    flood_timestamps = timestamps(~isEbb);

    F_struct_total = exceedance_probability(struct('Discharge', velocities, 'time', timestamps));
    F = F_struct_total.F;

    F_struct_ebb = exceedance_probability(struct('Discharge', ebb_velocities, 'time', ebb_timestamps));
    F_ebb = F_struct_ebb.F;

    F_struct_flood = exceedance_probability(struct('Discharge', flood_velocities, 'time', flood_timestamps));
    F_flood = F_struct_flood.F;

    bin_size = options.bin_size;

    % Remove duplicate points before interpolation
    [ebb_velocities, unique_idx_ebb] = unique(ebb_velocities);
    [flood_velocities, unique_idx_flood] = unique(flood_velocities);
    [velocities, unique_idx_velocities] = unique(velocities);

    F_ebb = F_ebb(unique_idx_ebb);
    F_flood = F_flood(unique_idx_flood);
    F = F(unique_idx_velocities);

    decimals = round(bin_size/0.1);
    s_new = round(min(velocities), decimals):bin_size:round(max(velocities), decimals)+bin_size;

    f_total = griddedInterpolant(velocities, F, 'linear', 'nearest');
    f_ebb = griddedInterpolant(ebb_velocities, F_ebb, 'linear', 'nearest');
    f_flood = griddedInterpolant(flood_velocities, F_flood, 'linear', 'nearest');

    F_total = f_total(s_new);
    F_ebb = f_ebb(s_new);
    F_flood = f_flood(s_new);

    % This differs from the python version but is necessary to keep the exceedance probability between 0 and 100
    F_max_total = max(max(F_ebb), max(F_flood)) * 2;

    area(s_new, [F_ebb/F_max_total*100; F_flood/F_max_total*100].', 'EdgeColor', 'none');

    xlabel('Velocity [\itm/s\rm]','FontSize',20);
    ylabel('Probability of Exceedance','FontSize',18);
    legend('Ebb','Flood')
    grid on
    title(options.title)

    len = strlength(options.savepath);
    if len > 1
        exportgraphics(gca, options.savepath);
    end

end

function isEbb = flood_or_ebb(directions, flood, ebb)
    max_angle = max(ebb, flood);
    min_angle = min(ebb, flood);

    lower_split = rem((min_angle + (360 - max_angle + min_angle)/2 ), 360);
    upper_split = lower_split + 180;

    if (lower_split <= ebb) && (ebb < upper_split)
        isEbb = ((directions < upper_split) & (directions >= lower_split));
    else
        isEbb = ~((directions < upper_split) & (directions >= lower_split));
    end
end
