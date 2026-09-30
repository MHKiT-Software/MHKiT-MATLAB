function M = mhkit_binned_statistic_2d(x, y, values, statistic, bin_spec, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute a 2D binned statistic on a bin grid and return it as a matrix struct
%
% Native MATLAB equivalent of scipy.stats.binned_statistic_2d.
%
% Bins are half-open [edge_i, edge_i+1) except the last bin, which is
% closed [edge_n-1, edge_n]. x, y, and values are paired by index:
% element k of each describes one sample. Samples with x or y outside
% the edges are dropped with a warning, but still count toward the
% 'probability' denominator. Samples with NaN in x, y, or values are
% dropped only if omitnan is true; otherwise a warning is issued and NaN
% values propagate into a bin.
%
% Parameters
% ------------
% x : vector
%   x position of each sample
% y : vector
%   y position of each sample, same length as x
% values : vector
%   Value of each sample to compute statistics on, same length as x
% statistic : char or string
%   The statistic to compute. One of 'mean', 'std', 'median', 'count',
%   'sum', 'min', 'max', 'probability', or 'frequency'. 'probability'
%   and 'frequency' are the same statistic.
%   'std' is the population standard deviation (1/N), following
%   MHKiT-Python convention.
%   'probability' is the count in each bin divided by the total number
%   of elements of values, including any samples outside the edges.
% bin_spec : struct
%   Bin spec for each axis, passed to mhkit_define_bins_2d. Each axis
%   struct must have exactly one of these three field sets.
%     bin_spec.x : struct
%       x axis spec, one of:
%         bin_spec.x.start : double
%           Lower edge of the first bin
%         bin_spec.x.stop : double
%           Upper edge of the last bin, stop - start a whole number of widths
%         bin_spec.x.width : double
%           Bin width, positive. Uniform bins from start to stop.
%       or
%         bin_spec.x.edges : numeric vector
%           Bin edges, strictly increasing, at least two elements.
%           Centers are the midpoints. Bins may be non-uniform.
%       or
%         bin_spec.x.centers : numeric vector
%           Bin centers, strictly increasing, at least two elements.
%           Edges are the midpoints between centers, extended half a
%           spacing beyond the end centers. Bins may be non-uniform.
%     bin_spec.y : struct
%       y axis spec, same three forms as bin_spec.x
% function_name : char (optional)
%   Name of the public function the user called, used in error and
%   warning identifiers so they read as coming from that function.
%   Default is this function's name.
% omitnan : logical (optional) default false
%   If true, a NaN at index k in any of x, y, or values removes element
%   k from all three together, so the remaining elements stay paired by
%   index and no element of values is ever binned by another index's x
%   or y. The removed elements do not count toward any bin or toward the
%   'probability' denominator. For example, if x(7) is NaN then y(7) and
%   values(7) are also removed, even though they are not NaN.
%   If false, a warning is issued when any NaN is present: a NaN element
%   of x or y falls in no bin but still counts toward the 'probability'
%   denominator, and a bin containing a NaN element of values returns
%   NaN for every statistic except 'count' and 'probability'.
%
% Returns
% ---------
% M : struct
%   M.values : matrix
%     Statistic per bin, y bins down the rows and x bins across the
%     columns, so the matrix reads like the plotted image. Empty bins
%     follow scipy: 0 for 'count', 'sum', and 'probability'; NaN
%     otherwise.
%   M.stat : char
%     Statistic used
%   M.x_bins : row vector
%     x bin centers
%   M.y_bins : row vector
%     y bin centers
%   M.x_edges : row vector
%     x bin edges, one more than the number of columns
%   M.y_edges : row vector
%     y bin edges, one more than the number of rows
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    x (:,1) {mustBeNumeric}
    y (:,1) {mustBeNumeric}
    values (:,1) {mustBeNumeric}
    statistic {mustBeTextScalar, mustBeMember(statistic, {'mean', 'std', 'median', 'count', 'sum', 'min', 'max', 'probability', 'frequency'})}
    bin_spec (1,1) struct
    options.function_name {mustBeTextScalar} = mfilename
    options.omitnan (1,1) logical = false
end

arguments (Output)
    M (1,1) struct
end

if length(x) ~= length(y) || length(x) ~= length(values)
    error(sprintf('MHKiT:%s:InvalidInput', options.function_name), ...
        'x, y, and values must have the same length');
end
if ~isfield(bin_spec, 'x') || ~isfield(bin_spec, 'y')
    error(sprintf('MHKiT:%s:InvalidInput', options.function_name), ...
        'bin_spec must be a struct with an x and a y bin spec, see mhkit_define_bins_2d');
end
grid = mhkit_define_bins_2d(bin_spec.x, bin_spec.y, 'function_name', options.function_name);
x_edges = grid.x.edges;
y_edges = grid.y.edges;

% One logical index over all three inputs, so a NaN at index k in any input
% removes element k from x, y, and values together and they stay aligned.
nan_index = isnan(x) | isnan(y) | isnan(values);
if any(nan_index)
    if options.omitnan
        x = x(~nan_index);
        y = y(~nan_index);
        values = values(~nan_index);
    else
        warning(sprintf('MHKiT:%s:NaNInput', options.function_name), ...
            ['%d of %d elements have NaN in x, y, or values. A NaN in x or y ' ...
             'falls in no bin and a NaN in values makes its bin NaN. Set ' ...
             'omitnan=true to remove these elements from all three inputs.'], ...
            nnz(nan_index), numel(nan_index));
    end
end

total_count = length(values);

outside = x < x_edges(1) | x > x_edges(end) | y < y_edges(1) | y > y_edges(end);
if any(outside)
    warning(sprintf('MHKiT:%s:OutsideEdges', options.function_name), ...
        ['%d of %d elements have x or y outside x_edges [%g, %g] or y_edges [%g, %g] ' ...
         'and are dropped. They still count toward the probability denominator.'], ...
        nnz(outside), numel(outside), x_edges(1), x_edges(end), y_edges(1), y_edges(end));
end

% Map the statistic name to a function and its empty-bin value
switch char(statistic)
    case 'mean'
        stat_func = @mean;
        null_value = NaN;
    case 'std'
        % Weight 1 gives the population standard deviation (1/N). MATLAB's
        % default is 1/(N-1), but MHKiT-Python convention is using scipy.stats.binned_statistic_2d 
        % which uses 1/N.
        stat_func = @(v) std(v, 1);
        null_value = NaN;
    case 'median'
        stat_func = @median;
        null_value = NaN;
    case 'count'
        stat_func = @numel;
        null_value = 0;
    case 'sum'
        stat_func = @sum;
        null_value = 0;
    case 'min'
        stat_func = @min;
        null_value = NaN;
    case 'max'
        stat_func = @max;
        null_value = NaN;
    case {'probability', 'frequency'}
        stat_func = @(v) numel(v) / total_count;
        null_value = 0;
end

nx = length(x_edges) - 1;
ny = length(y_edges) - 1;

% discretize returns NaN for elements outside the edges. Drop those so
% accumarray only sees valid subscripts, then apply the statistic per bin
% with null_value filling the empty bins. y first so y bins run down the
% rows and x bins across the columns.
x_index = discretize(x, x_edges);
y_index = discretize(y, y_edges);
in_bin = ~isnan(x_index) & ~isnan(y_index);
result = accumarray([y_index(in_bin), x_index(in_bin)], values(in_bin), [ny, nx], stat_func, null_value);

M = struct();
M.values = result;
M.stat = char(statistic);
M.x_bins = grid.x.centers';
M.y_bins = grid.y.centers';
M.x_edges = x_edges';
M.y_edges = y_edges';

end

