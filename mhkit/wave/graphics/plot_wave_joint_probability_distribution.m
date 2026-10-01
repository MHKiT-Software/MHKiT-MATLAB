function [ax, jpd] = plot_wave_joint_probability_distribution(Hm0, Te, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the joint probability distribution (JPD) of sea states, the
% fraction of records in each Hm0 and Te bin, as a percentage
%
% By default the bins follow IEC TS 62600-100: 0.5 m wide for Hm0 and
% 1 s wide for Te with edges on 0, 0.5, 1, ... m and 0, 1, 2, ... s,
% extending to the largest sea state. The returned matrix struct has
% the same layout as capture_width_matrix so its bins can be reused.
%
% Parameters
% ------------
%     Hm0 : vector
%         Significant wave height of each record [m]
%     Te : vector
%         Energy period of each record [s], same length as Hm0
%     Hm0_bins : vector (optional)
%         Hm0 bin centers [m]. Default IEC 0.5 m bins
%     Te_bins : vector (optional)
%         Te bin centers [s]. Default IEC 1 s bins
%     xlabel : string (optional)
%         Default "Energy Period, T_e [sec]"
%     ylabel : string (optional)
%         Default "Significant Wave Height, H_{m0} [m]"
%     title : string (optional)
%         Default "Joint Probability Distribution"
%     colormap : string or matrix (optional)
%         Colormap name ("viridis", any cmocean or MATLAB colormap name,
%         "-<colormap name>" or "<colormap name>_r" to flip) or an N x 3 RGB matrix, see mhkit_colormap.
%         Default "viridis"
%     annotate : logical (optional)
%         Print the percentage in each occupied bin. Default true
%     font_size : double (optional)
%         Font size of the bin labels in points. Default 8
%     trim_to_data : logical (optional)
%         Limit the axes to the occupied bins plus one empty bin around
%         them. Default true
%     ax : axes handle (optional)
%         Axes to plot into. Default is the current axes
%     savepath : string (optional)
%         Path and filename to save the figure. Default none
%
% Returns
% ---------
%     ax : axes handle
%         Axes the JPD was plotted into
%     jpd : structure
%         jpd.values : matrix
%             Fraction of records in each bin, Hm0 rows by Te columns,
%             summing to one
%         jpd.stat : char
%             'probability'
%         jpd.x_bins, x_edges : vector
%             Te bin centers and edges [s]
%         jpd.y_bins, y_edges : vector
%             Hm0 bin centers and edges [m]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric, mustBeVector}
    Te {mustBeNumeric, mustBeVector}
    options.Hm0_bins {mustBeNumeric} = []
    options.Te_bins {mustBeNumeric} = []
    options.xlabel {mustBeTextScalar} = "Energy Period, T_e [sec]"
    options.ylabel {mustBeTextScalar} = "Significant Wave Height, H_{m0} [m]"
    options.title {mustBeTextScalar} = "Joint Probability Distribution"
    options.colormap = "viridis"
    options.annotate (1,1) logical = true
    options.font_size (1,1) {mustBeNumeric, mustBePositive} = 8
    options.trim_to_data (1,1) logical = true
    options.ax = []
    options.savepath {mustBeTextScalar} = ""
end

arguments (Output)
    ax matlab.graphics.axis.Axes
    jpd struct
end

if numel(Hm0) ~= numel(Te)
    error('MHKiT:plot_wave_joint_probability_distribution:InvalidInput', ...
        'Hm0 and Te must have the same number of elements');
end

% IEC TS 62600-100 bin centers, so the edges fall on multiples of the width
if isempty(options.Hm0_bins)
    options.Hm0_bins = 0.25:0.5:ceil(max(Hm0) / 0.5) * 0.5;
end
if isempty(options.Te_bins)
    options.Te_bins = 0.5:1:ceil(max(Te));
end

bin_spec = struct('x', struct('centers', options.Te_bins), ...
                  'y', struct('centers', options.Hm0_bins));
% 'probability' ignores the values argument, so Hm0 is passed as a placeholder
jpd = mhkit_binned_statistic_2d(Te, Hm0, Hm0, 'probability', bin_spec, ...
    'function_name', mfilename);

% Plot as a percentage with the empty bins left blank
jpd_percent = jpd;
jpd_percent.values = 100 * jpd.values;
jpd_percent.values(jpd_percent.values == 0) = NaN;

if isempty(options.ax)
    options.ax = gca;
end
ax = mhkit_plot_matrix(jpd_percent, 'xlabel', options.xlabel, 'ylabel', options.ylabel, ...
    'zlabel', 'Occurrence [%]', 'colormap', options.colormap, ...
    'show_values', options.annotate, 'value_format', '%.2f %%', ...
    'font_size', options.font_size, 'trim_to_data', options.trim_to_data, ...
    'ax', options.ax);
title(ax, options.title);

if strlength(options.savepath) > 0
    saveas(ancestor(ax, 'figure'), options.savepath);
end

end
