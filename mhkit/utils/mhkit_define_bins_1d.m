function bins = mhkit_define_bins_1d(spec, function_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute bin centers and edges for one axis
%
% Parameters
% ------------
% spec : struct
%   Must have exactly one of the following three field sets.
%   Uniform bins beginning at start, each width wide:
%     spec.start : double
%       Lower edge of the first bin
%     spec.stop : double
%       Upper edge of the last bin. stop - start must be a whole number
%       of widths, so for start 0 and width 0.5, stop 3.5 is valid and
%       stop 3.7 is an error. Round a data maximum up to the width first.
%     spec.width : double
%       Bin width, positive
%   User-supplied edges, centers are the midpoints:
%     spec.edges : vector
%       Bin edges, strictly increasing, at least two elements
%   User-supplied centers, interior edges are the midpoints between
%   neighbors and the first and last edges extend half the neighboring
%   spacing beyond the end centers:
%     spec.centers : vector
%       Bin centers, strictly increasing, at least two elements
%
% function_name : char
%   Name of the public function the user called, used in error
%   identifiers so errors read as coming from that function, following
%   mhkit_standardize_spectrum_input
%
% Returns
% ---------
% bins : struct
%   bins.centers : column vector
%     Bin centers, one per bin
%   bins.edges : column vector
%     Bin edges, one more than the number of bins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spec (1,1) struct
    function_name {mustBeTextScalar}
end

arguments (Output)
    bins (1,1) struct
end

% Sort so the user can list the fields in any order, and make it a row
% cell so it can be compared to the literal field sets below
fields = sort(fieldnames(spec))';
if isequal(fields, {'start', 'stop', 'width'})
    for f = fields
        v = spec.(f{1});
        if ~(isnumeric(v) && isscalar(v) && isfinite(v))
            error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
                '%s must be a finite numeric scalar', f{1});
        end
    end
    if spec.width <= 0
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            'width (%g) must be positive', spec.width);
    end
    if spec.stop <= spec.start
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            'stop (%g) must be greater than start (%g)', spec.stop, spec.start);
    end
    n_widths = (spec.stop - spec.start) / spec.width;
    n_bins = round(n_widths);
    % Small tolerance because 0.3 / 0.1 is 3.0000000000000004 in floating point
    if abs(n_widths - n_bins) > 1e-9
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            ['stop - start (%g) must be a whole number of widths (%g). ' ...
             'Round stop up to the next multiple of width, or pass explicit ' ...
             'edges for a non-uniform grid.'], ...
            spec.stop - spec.start, spec.width);
    end
    % Transpose so edges and centers are column vectors
    edges = spec.start + spec.width * (0:n_bins)';
    centers = edges(1:end-1) + spec.width / 2;
elseif isequal(fields, {'edges'}) || isequal(fields, {'centers'})
    v = spec.(fields{1});
    if ~(isnumeric(v) && isvector(v) && numel(v) >= 2)
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s must be a numeric vector with at least two elements', fields{1});
    end
    % Force a column vector
    v = v(:);
    % Written as "not greater than" so a NaN, which fails every comparison,
    % is rejected along with repeated or decreasing values
    if any(~(diff(v) > 0))
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s must be strictly increasing with no NaN', fields{1});
    end
    if isequal(fields, {'edges'})
        edges = v;
        % Each center is the midpoint of its two edges
        centers = (edges(1:end-1) + edges(2:end)) / 2;
    else
        centers = v;
        % Interior edges are midpoints between neighboring centers
        mid = (centers(1:end-1) + centers(2:end)) / 2;
        % End edges extend half the neighboring spacing beyond the end centers
        first = centers(1) - (centers(2) - centers(1)) / 2;
        last = centers(end) + (centers(end) - centers(end-1)) / 2;
        edges = [first; mid; last];
    end
else
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        'spec must have exactly the fields start, stop, and width; or edges; or centers');
end

bins = struct();
bins.centers = centers;
bins.edges = edges;

end
