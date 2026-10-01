function WEFM = wave_energy_flux_matrix(Hm0, Te, J, statistic, Hm0_bins, Te_bins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a wave energy flux matrix for a given statistic
%
% Parameters
% ------------
% Hm0 : vector [m]
%   Significant wave height from spectra
% Te : vector [s]
%   Energy period from spectra
% J : vector [W/m]
%   Wave energy flux from spectra
% statistic : char or string
%   Statistic for each bin. Options: 'mean', 'std', 'median',
%   'count', 'sum', 'min', 'max', 'probability', or 'frequency'.
%   'probability' and 'frequency' are the same statistic.
%   'std' is the population standard deviation (1/N), following
%   MHKiT-Python convention.
% Hm0_bins : numeric vector [m]
%   Hm0 bin centers, strictly increasing. Edges are the midpoints between
%   centers, extended half a spacing beyond the end centers, following
%   MHKiT-Python convention. Bins may be non-uniform.
% Te_bins : numeric vector [s]
%   Te bin centers, same rules as Hm0_bins
%
% Returns
% ---------
% WEFM : struct
%   WEFM.values : matrix
%     Wave energy flux matrix, Hm0 bins down the rows and Te bins across the columns
%   WEFM.stat : string
%     Statistic used
%   WEFM.x_bins : row vector [s]
%     Te bin centers
%   WEFM.y_bins : row vector [m]
%     Hm0 bin centers
%   WEFM.x_edges : row vector [s]
%     Te bin edges, one more than the number of columns
%   WEFM.y_edges : row vector [m]
%     Hm0 bin edges, one more than the number of rows
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Hm0 {mustBeNumeric}
    Te {mustBeNumeric}
    J {mustBeNumeric}
    statistic {mustBeTextScalar, mustBeMember(statistic, {'mean', 'std', 'median', 'count', 'sum', 'min', 'max', 'probability', 'frequency'})}
    Hm0_bins {mustBeNumeric, mustBeVector}
    Te_bins {mustBeNumeric, mustBeVector}
end

arguments (Output)
    WEFM struct
end

% Te is x and Hm0 is y, following the IEC TS 62600-100 scatter diagram
bin_spec = struct('x', struct('centers', Te_bins), 'y', struct('centers', Hm0_bins));

WEFM = mhkit_binned_statistic_2d(Te, Hm0, J, statistic, bin_spec, 'function_name', mfilename);

if any(diff(WEFM.x_edges) > 1.0)
    warning('MHKiT:wave_energy_flux_matrix:BinSpacing', ...
        'Energy period bins are greater than the IEC TS 62600-100 limit of 1.0 seconds.');
end
if any(diff(WEFM.y_edges) > 0.5)
    warning('MHKiT:wave_energy_flux_matrix:BinSpacing', ...
        'Significant wave height bins are greater than the IEC TS 62600-100 limit of 0.5 meters.');
end

end
