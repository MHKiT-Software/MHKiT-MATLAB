function clm = capture_length_matrix(Hm0, Te, L, statistic, Hm0_bins, Te_bins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Deprecated alias for capture_width_matrix.
%
% IEC TS 62600-100 Ed. 2.0 replaces "capture length" with "capture
% width". This function will be removed in MHKiT-MATLAB v1.3; use
% capture_width_matrix instead.
%
% Parameters
% ------------
% Hm0 : vector [m]
%   Significant wave height from spectra
% Te : vector [s]
%   Energy period from spectra
% L : vector [m]
%   Capture length
% statistic : char or string
%   Statistic for each bin. Options: 'mean', 'std', 'median',
%   'count', 'sum', 'min', 'max', 'probability', or 'frequency'.
%   'probability' and 'frequency' are the same statistic.
% Hm0_bins : numeric vector [m]
%   Hm0 bin centers, see capture_width_matrix
% Te_bins : numeric vector [s]
%   Te bin centers, see capture_width_matrix
%
% Returns
% ---------
% clm : struct
%   Same layout as capture_width_matrix: values (Hm0 rows x Te columns),
%   stat, x_bins, y_bins, x_edges, y_edges with Te on x and Hm0 on y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric}
    Te {mustBeNumeric}
    L {mustBeNumeric}
    statistic {mustBeTextScalar}
    Hm0_bins {mustBeNumeric, mustBeVector}
    Te_bins {mustBeNumeric, mustBeVector}
end

arguments (Output)
    clm struct
end

warning('MHKiT:capture_length_matrix:DeprecatedFunction', ...
    ['IEC TS 62600-100 Ed. 2.0 replaces "capture length" with "capture ' ...
    'width". capture_length_matrix will be removed in MHKiT-MATLAB v1.3. ' ...
    'Use capture_width_matrix instead.']);

clm = capture_width_matrix(Hm0, Te, L, statistic, Hm0_bins, Te_bins);

end
