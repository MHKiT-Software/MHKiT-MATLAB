function ax = plot_matrix(M, Mtype, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots a wave performance matrix with Hm0 on the y axis and Te on the x axis
%
% Parameters
% ------------
% M : struct
%   Matrix from capture_width_matrix, wave_energy_flux_matrix, or power_matrix
%     M.values : matrix
%       One value per bin, Hm0 bins down the rows and Te bins across the columns
%     M.stat : string
%       Statistic used, shown in the title
%     M.x_edges : vector [s]
%       Te bin edges
%     M.y_edges : vector [m]
%       Hm0 bin edges
% Mtype : string
%   Type of matrix (e.g. "Capture Width", "Power") used in the plot title
% savepath : string (optional)
%   Path and filename to save the figure. Default none
% annotate : logical (optional)
%   Print each bin value in its cell. Default true
% value_format : string (optional)
%   sprintf format for the bin labels, e.g. '%.2f m'. Default '%.2f'
% font_size : double (optional)
%   Font size of the bin labels in points. Default 7
% trim_to_data : logical (optional)
%   Limit the axes to the bins with data plus one empty bin around them. Default false
% zlabel : string (optional)
%   Colorbar label, e.g. "Capture Width [m]". Default none
%
% Returns
% ---------
% ax : matlab.graphics.axis.Axes
%   Axes containing the plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    M (1,1) struct
    Mtype {mustBeTextScalar}
    options.savepath {mustBeTextScalar} = ""
    options.annotate (1,1) logical = true
    options.zlabel {mustBeTextScalar} = ""
    options.trim_to_data (1,1) logical = false
    options.value_format {mustBeTextScalar} = '%.2f'
    options.font_size (1,1) {mustBeNumeric, mustBePositive} = 7
end

arguments (Output)
    ax (1,1) matlab.graphics.axis.Axes
end

ax = mhkit_plot_matrix(M, 'xlabel', 'Energy Period, T_e [sec]', 'ylabel', 'Significant Wave Height, H_{m0} [m]', ...
    'show_values', options.annotate, 'zlabel', options.zlabel, ...
    'trim_to_data', options.trim_to_data, 'value_format', options.value_format, ...
    'font_size', options.font_size);
title(ax, string(Mtype) + ": " + string(M.stat));

if strlength(options.savepath) > 0
    saveas(ancestor(ax, 'figure'), options.savepath);
end

end
