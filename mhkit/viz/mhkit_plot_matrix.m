function ax = mhkit_plot_matrix(M, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot a binned matrix visualization, one colored cell per bin
%
% Cells are drawn between the bin edges so every bin is visible, with the
% tick marks at the bin centers and empty (NaN) bins left blank.
%
% Parameters
% ------------
% M : struct
%   Matrix struct from mhkit_binned_statistic_2d, capture_width_matrix,
%   wave_energy_flux_matrix, or power_matrix
%     M.values : matrix
%       One value per bin, y bins down the rows and x bins across the columns
%     M.x_edges : vector
%       x bin edges, one more than the number of columns
%     M.y_edges : vector
%       y bin edges, one more than the number of rows
% xlabel : string (optional)
%   x axis label. Default "Te"
% ylabel : string (optional)
%   y axis label. Default "Hm0"
% zlabel : string (optional)
%   Colorbar label. Default none
% show_values : logical (optional)
%   Print each bin value to two decimals in its cell. Default true
% ax : axes handle (optional)
%   Axes to plot into. Default is a new figure
% savepath : string (optional)
%   Path and filename to save the figure. Default none
%
% Returns
% ---------
% ax : matlab.graphics.axis.Axes
%   Axes containing the plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    M (1,1) struct
    options.xlabel {mustBeTextScalar} = "Te"
    options.ylabel {mustBeTextScalar} = "Hm0"
    options.zlabel {mustBeTextScalar} = ""
    options.show_values (1,1) logical = true
    options.ax = []
    options.savepath {mustBeTextScalar} = ""
end

arguments (Output)
    ax (1,1) matlab.graphics.axis.Axes
end

for f = {'values', 'x_edges', 'y_edges'}
    if ~isfield(M, f{1})
        error('MHKiT:mhkit_plot_matrix:InvalidInput', ...
            'M must have values, x_edges, and y_edges fields, missing %s', f{1});
    end
end
x_edges = M.x_edges(:)';
y_edges = M.y_edges(:)';
[rows, cols] = size(M.values);
if rows ~= numel(y_edges) - 1 || cols ~= numel(x_edges) - 1
    error('MHKiT:mhkit_plot_matrix:InvalidInput', ...
        'values is %dx%d but the edges define %d rows and %d columns', ...
        rows, cols, numel(y_edges) - 1, numel(x_edges) - 1);
end

if isempty(options.ax)
    figure;
    ax = gca;
elseif isa(options.ax, 'matlab.graphics.axis.Axes')
    ax = options.ax;
else
    error('MHKiT:mhkit_plot_matrix:InvalidInput', 'ax must be an axes handle');
end

% pcolor colors cell (i,j) with C(i,j) and ignores the last row and column
% of C, so pad values by one so every bin is drawn between its edges
padded = [M.values, nan(rows, 1); nan(1, cols + 1)];
pcolor(ax, x_edges, y_edges, padded);
shading(ax, 'flat');
colormap(ax, viridis_colormap());
cb = colorbar(ax);
if strlength(options.zlabel) > 0
    cb.Label.String = options.zlabel;
end

x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
xticks(ax, x_centers);
yticks(ax, y_centers);
xlim(ax, x_edges([1 end]));
ylim(ax, y_edges([1 end]));
xlabel(ax, options.xlabel);
ylabel(ax, options.ylabel);

if options.show_values
    % Dark text on the bright upper half of the colormap, white on the rest
    mid = mean(clim(ax));
    for i = 1:rows
        for j = 1:cols
            v = M.values(i, j);
            if ~isnan(v)
                if v > mid
                    color = [0 0 0];
                else
                    color = [1 1 1];
                end
                text(ax, x_centers(j), y_centers(i), sprintf('%.2f', v), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'Color', color);
            end
        end
    end
end

if strlength(options.savepath) > 0
    saveas(ancestor(ax, 'figure'), options.savepath);
end

end
