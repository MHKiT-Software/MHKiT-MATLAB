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
end

arguments (Output)
    ax (1,1) matlab.graphics.axis.Axes
end

ax = mhkit_plot_matrix(M, 'xlabel', 'Te [s]', 'ylabel', 'Hm0 [m]', ...
    'show_values', options.annotate);
title(ax, string(Mtype) + ": " + string(M.stat));

if strlength(options.savepath) > 0
    saveas(ancestor(ax, 'figure'), options.savepath);
end

end
