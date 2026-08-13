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
% statistic : string or function_handle
%   Statistic for each bin. Options: 'mean', 'std', 'median',
%   'count', 'sum', 'min', 'max', 'frequency', or a function handle.
% Hm0_bins : vector [m]
%   Bin centers for Hm0
% Te_bins : vector [s]
%   Bin centers for Te
%
% Returns
% ---------
% clm : struct
%   clm.values : matrix
%     Capture length matrix (Hm0_bins x Te_bins)
%   clm.stat : string
%     Statistic used
%   clm.Hm0_bins : vector [m]
%     Hm0 bin centers
%   clm.Te_bins : vector [s]
%     Te bin centers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric}
    Te {mustBeNumeric}
    L {mustBeNumeric}
    statistic
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
