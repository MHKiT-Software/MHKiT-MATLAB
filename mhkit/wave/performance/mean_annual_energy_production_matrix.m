function maep = mean_annual_energy_production_matrix(CWM, JM, frequency)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates mean annual energy production (MAEP) from matrix data
%
% MAEP = T * nansum(CWM * JM * frequency)
% where T = 8766 hours (average length of a year)
%
% Parameters
% ------------
% CWM : struct or matrix
%   Capture width matrix. If struct, uses CWM.values
% JM : struct or matrix
%   Wave energy flux matrix. If struct, uses JM.values
% frequency : struct or matrix
%   Data frequency for each bin (must sum to 1). If struct, uses frequency.values
%
% Returns
% ---------
% maep : double [W*h]
%   Mean annual energy production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    CWM
    JM
    frequency
end

arguments (Output)
    maep {mustBeNumeric}
end

% Extract values from structs if needed
if isstruct(CWM)
    CWM_vals = CWM.values;
else
    CWM_vals = CWM;
end

if isstruct(JM)
    JM_vals = JM.values;
else
    JM_vals = JM;
end

if isstruct(frequency)
    freq_vals = frequency.values;
else
    freq_vals = frequency;
end

% Validate dimensions match
if ~isequal(size(CWM_vals), size(JM_vals)) || ~isequal(size(CWM_vals), size(freq_vals))
    error('MHKiT:mean_annual_energy_production_matrix:DimensionMismatch', ...
        'CWM, JM, and frequency must have the same dimensions');
end

% Validate frequency sums to 1
freq_sum = sum(freq_vals(:), 'omitnan');
if abs(freq_sum - 1) > 1e-6
    error('MHKiT:mean_annual_energy_production_matrix:InvalidFrequency', ...
        'Frequency components must sum to one. Got: %g', freq_sum);
end

T = 8766;  % Average length of a year in hours

maep = T * sum(CWM_vals(:) .* JM_vals(:) .* freq_vals(:), 'omitnan');

end
