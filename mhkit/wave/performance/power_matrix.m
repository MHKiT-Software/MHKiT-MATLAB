function PM = power_matrix(CWM, JM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates a power matrix from a capture width matrix and wave energy flux matrix
%
% PM = CWM * JM (element-wise multiplication)
%
% Parameters
% ------------
% CWM : struct
%   Capture width matrix from capture_width_matrix:
%     CWM.values : matrix
%       Capture width values
%     CWM.stat : string
%       Statistic used (e.g., 'mean')
%     CWM.x_bins, CWM.y_bins : row vectors
%       Te and Hm0 bin centers
%     CWM.x_edges, CWM.y_edges : row vectors
%       Te and Hm0 bin edges, must match JM
% JM : struct
%   Wave energy flux matrix structure:
%     JM.values : matrix
%       Wave energy flux values
%     JM.x_bins, JM.y_bins : row vectors
%       Te and Hm0 bin centers
%     JM.x_edges, JM.y_edges : row vectors
%       Te and Hm0 bin edges, must match CWM
%
% Returns
% ---------
% PM : struct
%   Power matrix structure:
%     PM.values : matrix
%       Power matrix values
%     PM.stat : string
%       Statistic from CWM
%     PM.x_bins, PM.y_bins : row vectors
%       Te and Hm0 bin centers from CWM
%     PM.x_edges, PM.y_edges : row vectors
%       Te and Hm0 bin edges from CWM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    CWM struct
    JM struct
end

arguments (Output)
    PM struct
end

% Both inputs must be matrix structs from mhkit_binned_statistic_2d on the same grid
for f = {'values', 'stat', 'x_bins', 'y_bins', 'x_edges', 'y_edges'}
    if ~isfield(CWM, f{1}) || ~isfield(JM, f{1})
        error('MHKiT:power_matrix:InvalidInput', ...
            'CWM and JM must both be matrix structs from capture_width_matrix and wave_energy_flux_matrix, missing field %s', f{1});
    end
end
if ~isequal(CWM.x_edges, JM.x_edges) || ~isequal(CWM.y_edges, JM.y_edges)
    error('MHKiT:power_matrix:InvalidInput', ...
        'CWM and JM must be binned on the same x and y edges');
end

% Build output structure
PM = struct();
PM.values = CWM.values .* JM.values;
PM.stat = CWM.stat;
PM.x_bins = CWM.x_bins;
PM.y_bins = CWM.y_bins;
PM.x_edges = CWM.x_edges;
PM.y_edges = CWM.y_edges;

end
