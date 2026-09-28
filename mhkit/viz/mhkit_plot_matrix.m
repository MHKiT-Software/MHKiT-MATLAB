function ax = mhkit_plot_matrix(M, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot a binned matrix visualization, one colored cell per bin
%
% Cells are drawn between the bin edges so every bin is visible, with the
% tick marks and a light grid at the bin edges and empty (NaN) bins left
% blank.
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
% trim_to_data : logical (optional)
%   Limit the axes to the bins that hold data plus one empty bin on every
%   side, instead of the full grid. Default false
% show_values : logical (optional)
%   Print a label in each bin that holds data. Default true
% value_format : string (optional)
%   sprintf format for the bin value, e.g. '%.1f kW'. Default '%.2f'
% labels : cell array (optional)
%   Custom label per bin, same size as M.values, used instead of
%   value_format when given. Empty labels are skipped. Default none
% font_size : double (optional)
%   Font size of the bin labels in points. Default 7
% ax : axes handle (optional)
%   Axes to plot into. Default is the current axes, which opens a new
%   figure only if none exists, like the built-in plot functions
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
    options.trim_to_data (1,1) logical = false
    options.show_values (1,1) logical = true
    options.value_format {mustBeTextScalar} = '%.2f'
    options.labels cell = {}
    options.font_size (1,1) {mustBeNumeric, mustBePositive} = 7
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
    ax = gca;
elseif isa(options.ax, 'matlab.graphics.axis.Axes')
    ax = options.ax;
else
    error('MHKiT:mhkit_plot_matrix:InvalidInput', 'ax must be an axes handle');
end

% pcolor colors cell (i,j) with C(i,j) and ignores the last row and column
% of C, so pad values by one so every bin is drawn between its edges
padded = [M.values, nan(rows, 1); nan(1, cols + 1)];
h = pcolor(ax, x_edges, y_edges, padded);
% Light grid on every cell edge so the bins read as bins, including empty ones
h.EdgeColor = [0.8 0.8 0.8];
h.LineWidth = 0.25;
colormap(ax, viridis_colormap());
cb = colorbar(ax);
if strlength(options.zlabel) > 0
    cb.Label.String = options.zlabel;
end

x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
xticks(ax, x_edges);
yticks(ax, y_edges);
if options.trim_to_data && any(~isnan(M.values(:)))
    % Bounding box of the occupied bins, padded by one bin and clamped to the grid
    occupied_cols = find(any(~isnan(M.values), 1));
    occupied_rows = find(any(~isnan(M.values), 2));
    xlim(ax, x_edges([max(occupied_cols(1) - 1, 1), min(occupied_cols(end) + 2, cols + 1)]));
    ylim(ax, y_edges([max(occupied_rows(1) - 1, 1), min(occupied_rows(end) + 2, rows + 1)]));
else
    xlim(ax, x_edges([1 end]));
    ylim(ax, y_edges([1 end]));
end
xlabel(ax, options.xlabel);
ylabel(ax, options.ylabel);

if options.show_values
    use_labels = ~isempty(options.labels);
    if use_labels && ~isequal(size(options.labels), size(M.values))
        error('MHKiT:mhkit_plot_matrix:InvalidInput', ...
            'labels must be a cell array the same size as M.values');
    end
    % Dark text on the bright upper half of the colormap, white on the rest
    mid = mean(clim(ax));
    for i = 1:rows
        for j = 1:cols
            v = M.values(i, j);
            if isnan(v)
                continue
            end
            if use_labels
                label = options.labels{i, j};
            else
                % Round small nonzero values up so they never print as zero
                if v ~= 0 && abs(v) < 0.005
                    v = sign(v) * 0.01;
                end
                label = sprintf(options.value_format, v);
            end
            if isempty(label)
                continue
            end
            if M.values(i, j) > mid
                color = [0 0 0];
            else
                color = [1 1 1];
            end
            text(ax, x_centers(j), y_centers(i), label, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', color, 'FontSize', options.font_size);
        end
    end
end

if strlength(options.savepath) > 0
    saveas(ancestor(ax, 'figure'), options.savepath);
end

end
