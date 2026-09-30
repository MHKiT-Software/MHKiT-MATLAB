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
%   Limit the axes to the bins that hold data plus padding empty bins on
%   every side, instead of the full grid. Default false
% padding : integer (optional)
%   Number of empty bins shown around the occupied bins when trim_to_data
%   is true. Where the grid ends before the padding does, it is extended
%   with empty bins of the same width, except below an edge at zero since
%   the binned quantities are normally nonnegative. Default 1
% show_values : logical (optional)
%   Print a label in each bin that holds data. Default true
% value_format : string (optional)
%   sprintf format for the bin value, e.g. '%.1f kW'. Default '%.2f'
% labels : cell array (optional)
%   Custom label per bin, same size as M.values, used instead of
%   value_format when given. Empty labels are skipped. Default none
% font_size : double (optional)
%   Font size of the bin labels in points. Default 7
% colormap : string or matrix (optional)
%   Colormap name ("viridis", any cmocean or MATLAB colormap name,
%   "-<colormap name>" or "<colormap name>_r" to flip) or an N x 3 RGB matrix, see mhkit_colormap. Default "viridis"
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
    options.padding (1,1) {mustBeInteger, mustBeNonnegative} = 1
    options.show_values (1,1) logical = true
    options.value_format {mustBeTextScalar} = '%.2f'
    options.labels cell = {}
    options.font_size (1,1) {mustBeNumeric, mustBePositive} = 7
    options.colormap = "viridis"
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

values = M.values;
if options.trim_to_data && any(~isnan(values(:)))
    % Index range of the occupied bins plus the padding on every side
    occupied_cols = find(any(~isnan(values), 1));
    occupied_rows = find(any(~isnan(values), 2));
    first_col = occupied_cols(1) - options.padding;
    last_col = occupied_cols(end) + options.padding;
    first_row = occupied_rows(1) - options.padding;
    last_row = occupied_rows(end) + options.padding;
    % Extend the grid with empty bins where the padding runs past it
    [x_edges, values, first_col, last_col, col_shift] = extend_grid(x_edges, values, first_col, last_col, 2);
    [y_edges, values, first_row, last_row, row_shift] = extend_grid(y_edges, values, first_row, last_row, 1);
    x_limits = x_edges([first_col, last_col + 1]);
    y_limits = y_edges([first_row, last_row + 1]);
else
    x_limits = x_edges([1 end]);
    y_limits = y_edges([1 end]);
    row_shift = 0;
    col_shift = 0;
end
[rows, cols] = size(values);

% pcolor colors cell (i,j) with C(i,j) and ignores the last row and column
% of C, so pad values by one so every bin is drawn between its edges
padded = [values, nan(rows, 1); nan(1, cols + 1)];
h = pcolor(ax, x_edges, y_edges, padded);
% Light grid on every cell edge so the bins read as bins, including empty ones
h.EdgeColor = [0.8 0.8 0.8];
h.LineWidth = 0.25;
colormap(ax, mhkit_colormap(options.colormap));
cb = colorbar(ax);
if strlength(options.zlabel) > 0
    cb.Label.String = options.zlabel;
end

x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
xticks(ax, x_edges);
yticks(ax, y_edges);
% Keep tick labels upright; MATLAB rotates them automatically when they crowd
ax.XTickLabelRotation = 0;
ax.YTickLabelRotation = 0;
xlim(ax, x_limits);
ylim(ax, y_limits);
xlabel(ax, options.xlabel);
ylabel(ax, options.ylabel);

if options.show_values
    use_labels = ~isempty(options.labels);
    if use_labels && ~isequal(size(options.labels), size(M.values))
        error('MHKiT:mhkit_plot_matrix:InvalidInput', ...
            'labels must be a cell array the same size as M.values');
    end
    if use_labels
        % Place the labels on the extended grid
        labels = cell(size(values));
        labels(row_shift + (1:size(M.values, 1)), col_shift + (1:size(M.values, 2))) = options.labels;
    end
    % Dark text on the bright upper half of the colormap, white on the rest
    mid = mean(clim(ax));
    for i = 1:rows
        for j = 1:cols
            v = values(i, j);
            if isnan(v)
                continue
            end
            if use_labels
                label = labels{i, j};
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
            if values(i, j) > mid
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


function [edges, values, first, last, n_before] = extend_grid(edges, values, first, last, dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Extend one axis of the grid with empty bins so the bin index range
%    first:last exists, keeping the width of the end bins
%
%    The grid is not extended below an edge at zero, so first is clamped
%    there instead.
%
% Parameters
% ------------
%     edges : row vector
%         Bin edges along the axis
%     values : matrix
%         Bin values, y rows by x columns
%     first : integer
%         Wanted first bin index, may be less than 1
%     last : integer
%         Wanted last bin index, may exceed the number of bins
%     dim : integer
%         1 to extend the rows (y axis), 2 to extend the columns (x axis)
%
% Returns
% ---------
%     edges : row vector
%         Extended bin edges
%     values : matrix
%         Values padded with NaN bins
%     first : integer
%         First bin index in the extended grid
%     last : integer
%         Last bin index in the extended grid
%     n_before : integer
%         Number of empty bins added before the first original bin
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_bins = numel(edges) - 1;
n_before = max(1 - first, 0);
if edges(1) == 0
    n_before = 0;
    first = max(first, 1);
end
n_after = max(last - n_bins, 0);
if n_before > 0
    width = edges(2) - edges(1);
    edges = [edges(1) - width * (n_before:-1:1), edges];
end
if n_after > 0
    width = edges(end) - edges(end-1);
    edges = [edges, edges(end) + width * (1:n_after)];
end
pad_shape = size(values);
pad_shape(dim) = n_before;
before = nan(pad_shape);
pad_shape(dim) = n_after;
after = nan(pad_shape);
values = cat(dim, before, values, after);
first = first + n_before;
last = last + n_before;

end
