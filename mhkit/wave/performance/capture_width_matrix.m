function cwm = capture_width_matrix(Hm0, Te, CW, statistic, Hm0_bins, Te_bins)

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
% statistic : string or function_handle
%   Statistic for each bin. Options: 'mean', 'std', 'median',
%   'count', 'sum', 'min', 'max', 'frequency', or a function handle.
%   Note that 'std' uses a degree of freedom of N in accordance with
%   Formula D.5 of IEC TS 62600-100 Ed. 2.0 en 2024.
% Hm0_bins : vector [m]
%   Bin centers for Hm0
% Te_bins : vector [s]
%   Bin centers for Te
%
% Returns
% ---------
% cwm : struct
%   cwm.values : matrix
%     Capture width matrix (Hm0_bins x Te_bins)
%   cwm.stat : string
%     Statistic used
%   cwm.Hm0_bins : vector [m]
%     Hm0 bin centers
%   cwm.Te_bins : vector [s]
%     Te bin centers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric}
    Te {mustBeNumeric}
    CW {mustBeNumeric}
    statistic
    Hm0_bins {mustBeNumeric, mustBeVector}
    Te_bins {mustBeNumeric, mustBeVector}
end

arguments (Output)
    cwm struct
end

% Ensure column vectors
Hm0 = Hm0(:);
Te = Te(:);
CW = CW(:);
Hm0_bins = Hm0_bins(:);
Te_bins = Te_bins(:);

% Convert bin centers to edges
Hm0_edges = bin_centers_to_edges(Hm0_bins);
Te_edges = bin_centers_to_edges(Te_bins);

% Check bin spacing against IEC limits
Hm0_spacing = diff(Hm0_edges);
if any(Hm0_spacing > 0.5)
    warning('MHKiT:capture_width_matrix:BinSpacing', ...
        'Significant wave height bins are greater than the IEC TS 62600-100 limit of 0.5 meters.');
end

Te_spacing = diff(Te_edges);
if any(Te_spacing > 1.0)
    warning('MHKiT:capture_width_matrix:BinSpacing', ...
        'Energy period bins are greater than the IEC TS 62600-100 limit of 1.0 seconds.');
end

% Handle std with ddof=0 (population std) per IEC/TS 62600-100
if ischar(statistic) || isstring(statistic)
    if strcmpi(statistic, 'std')
        % Use population std (ddof=0) per IEC Formula D.5
        stat_func = @(x) std_ddof0(x);
    else
        stat_func = statistic;
    end
else
    stat_func = statistic;
end

% Compute binned statistic
values = binned_statistic_2d(Hm0, Te, CW, Hm0_edges, Te_edges, stat_func);

% Build output structure
cwm = struct();
cwm.values = values;
if ischar(statistic) || isstring(statistic)
    cwm.stat = char(statistic);
else
    cwm.stat = 'custom';
end
cwm.Hm0_bins = Hm0_bins';
cwm.Te_bins = Te_bins';

end

function s = std_ddof0(x)
    % Standard deviation with degree of freedom = N (population std)
    % Per IEC/TS 62600-100 Ed. 2.0 Formula D.5
    if isempty(x)
        s = NaN;
    elseif length(x) == 1
        s = 0;
    else
        s = std(x, 1);  % MATLAB std with w=1 gives population std
    end
end
