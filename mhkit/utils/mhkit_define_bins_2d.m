function bins = mhkit_define_bins_2d(x, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define a 2D bin grid from a spec struct along each axis
%
% Each axis is a struct with exactly one of three field sets, so the
% function knows which was supplied. Bins are computed by
% mhkit_define_bins_1d. For wave performance matrices, x is Te and y
% is Hm0, so x is the horizontal axis of the plotted matrix.
%
% Parameters
% ------------
% x : struct
%   Exactly one of the following three field sets.
%   Uniform bins beginning at start, each width wide, so the edges are
%   start, start + width, start + 2*width, ... and the centers sit
%   halfway between edges. For start 0, stop 2, width 0.5 the edges are
%   0, 0.5, 1, 1.5, 2 and the centers are 0.25, 0.75, 1.25, 1.75.
%     x.start : double
%       Lower edge of the first x bin
%     x.stop : double
%       Upper edge of the last x bin. stop - start must be a whole
%       number of widths. Round a data maximum up to the width first.
%     x.width : double
%       x bin width, positive.
%   User-supplied edges, centers are the midpoints:
%     x.edges : vector
%       x bin edges, strictly increasing, at least two elements
%   User-supplied centers, interior edges are the midpoints between
%   neighbors and the first and last edges extend half the neighboring
%   spacing beyond the end centers:
%     x.centers : vector
%       x bin centers, strictly increasing, at least two elements
% y : struct
%   Exactly one of the same three field sets, same rules as x.
%
% Returns
% ---------
% bins : struct
%   bins.x.centers : column vector
%     x bin centers, one per bin
%   bins.x.edges : column vector
%     x bin edges, one more than the number of bins
%   bins.y.centers : column vector
%     y bin centers
%   bins.y.edges : column vector
%     y bin edges
%
% Examples
% --------
%     >> % IEC TS 62600-100 typical grid is 1 s Te bins and 0.5 m Hm0 bins.
%     >> te_spec = struct('start', 2, 'stop', 15, 'width', 1);
%     >> hm0_spec = struct('start', 0, 'stop', 2, 'width', 0.5);
%     >> bins = mhkit_define_bins_2d(te_spec, hm0_spec);
%     >> bins.x.edges'
%     ans =
%          2     3     4     5     6     7     8     9    10    11    12    13    14    15
%
%     >> bins.y.edges'
%     ans =
%              0    0.5000    1.0000    1.5000    2.0000
%
%     >> % Non-uniform grid from explicit edges and centers
%     >> bins = mhkit_define_bins_2d(struct('edges', [2 4 8 15]), struct('centers', [0.25 0.75 1.25]));
%     >> bins.x.centers'
%     ans =
%          3     6    11.5000
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    x (1,1) struct
    y (1,1) struct
end

arguments (Output)
    bins (1,1) struct
end

bins = struct();
bins.x = mhkit_define_bins_1d(x, mfilename);
bins.y = mhkit_define_bins_1d(y, mfilename);

end
