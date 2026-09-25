function cwm = capture_width_matrix(Hm0, Te, CW, statistic, bin_spec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a capture width matrix for a given statistic
%
% Note that IEC TS 62600-100 Ed. 2.0 en 2024 section 9.2.4 requires
% capture width matrices for the mean, std, count, min, and max.
%
% Parameters
% ------------
% Hm0 : vector [m]
%   Significant wave height from spectra
% Te : vector [s]
%   Energy period from spectra
% CW : vector [m]
%   Capture width
% statistic : char or string
%   Statistic for each bin. Options: 'mean', 'std', 'median',
%   'count', 'sum', 'min', 'max', 'probability', or 'frequency'.
%   'probability' and 'frequency' are the same statistic.
%   'std' is the population standard deviation (1/N), following
%   MHKiT-Python convention.
% bin_spec : struct
%   Bin spec for each axis, see mhkit_define_bins_2d. Te is x and Hm0 is
%   y. Each axis struct must have exactly one of these three field sets.
%     bin_spec.x : struct
%       Te spec [s], one of:
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
%       Hm0 spec [m], same three forms as bin_spec.x
%
% Returns
% ---------
% cwm : struct
%   cwm.values : matrix
%     Capture width matrix, Hm0 bins down the rows and Te bins across the columns
%   cwm.stat : string
%     Statistic used
%   cwm.x_bins : row vector [s]
%     Te bin centers
%   cwm.y_bins : row vector [m]
%     Hm0 bin centers
%   cwm.x_edges : row vector [s]
%     Te bin edges, one more than the number of columns
%   cwm.y_edges : row vector [m]
%     Hm0 bin edges, one more than the number of rows
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric}
    Te {mustBeNumeric}
    CW {mustBeNumeric}
    statistic {mustBeTextScalar, mustBeMember(statistic, {'mean', 'std', 'median', 'count', 'sum', 'min', 'max', 'probability', 'frequency'})}
    bin_spec (1,1) struct
end

arguments (Output)
    cwm struct
end

cwm = mhkit_binned_statistic_2d(Te, Hm0, CW, statistic, bin_spec, 'function_name', mfilename);

% Te is x and Hm0 is y, following the IEC TS 62600-100 convention
if any(diff(cwm.x_edges) > 1.0)
    warning('MHKiT:capture_width_matrix:BinSpacing', ...
        'Energy period bins are greater than the IEC TS 62600-100 limit of 1.0 seconds.');
end
if any(diff(cwm.y_edges) > 0.5)
    warning('MHKiT:capture_width_matrix:BinSpacing', ...
        'Significant wave height bins are greater than the IEC TS 62600-100 limit of 0.5 meters.');
end

end
