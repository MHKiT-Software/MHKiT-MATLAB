function PM = power_matrix(LM, JM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a power matrix from a capture length matrix and wave energy flux matrix
%
% PM = LM * JM (element-wise multiplication)
%
% Parameters
% ------------
% LM : struct
%   Capture length matrix structure:
%     LM.values : matrix
%       Capture length values
%     LM.stat : string
%       Statistic used (e.g., 'mean')
%     LM.Hm0_bins : vector [m]
%       Hm0 bin centers
%     LM.Te_bins : vector [s]
%       Te bin centers
% JM : struct
%   Wave energy flux matrix structure:
%     JM.values : matrix
%       Wave energy flux values
%     JM.Hm0_bins : vector [m]
%       Hm0 bin centers
%     JM.Te_bins : vector [s]
%       Te bin centers
%
% Returns
% ---------
% PM : struct
%   Power matrix structure:
%     PM.values : matrix
%       Power matrix values
%     PM.stat : string
%       Statistic from LM
%     PM.Hm0_bins : vector [m]
%       Hm0 bin centers
%     PM.Te_bins : vector [s]
%       Te bin centers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    LM struct
    JM struct
end

arguments (Output)
    PM struct
end

% Validate input structures have values field
if ~isfield(LM, 'values')
    error('MHKiT:power_matrix:InvalidInput', ...
        'LM must be a structure with a values field');
end
if ~isfield(JM, 'values')
    error('MHKiT:power_matrix:InvalidInput', ...
        'JM must be a structure with a values field');
end

% Validate dimensions match
if ~isequal(size(LM.values), size(JM.values))
    error('MHKiT:power_matrix:DimensionMismatch', ...
        'LM.values and JM.values must have the same dimensions');
end

% Build output structure
PM = struct();
PM.values = LM.values .* JM.values;

if isfield(LM, 'stat')
    PM.stat = LM.stat;
else
    PM.stat = 'computed';
end

if isfield(LM, 'Hm0_bins')
    PM.Hm0_bins = LM.Hm0_bins;
end

if isfield(LM, 'Te_bins')
    PM.Te_bins = LM.Te_bins;
end

end
