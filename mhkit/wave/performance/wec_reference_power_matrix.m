function M = wec_reference_power_matrix(device)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Power matrix of a DOE Reference Model Wave Energy Converter
%
% The Reference Model Project, sponsored by the U.S. Department of Energy,
% developed open-source marine energy point designs as reference models to
% benchmark technology performance and costs:
% https://openei.org/wiki/PRIMRE/Signature_Projects/Reference_Model
% Returns the device power [kW] over significant wave height and energy
% period bins for Reference Model 3 (wave point absorber), 5 (oscillating
% surge flap), or 6 (oscillating water column), as distributed with the
% National Laboratory of the Rockies (NLR) System Advisor Model (SAM,
% BSD-3-Clause, https://github.com/NatLabRockies/SAM). The struct has the same layout
% as the matrices from capture_width_matrix and power_matrix, so it can be
% plotted with mhkit_plot_matrix and used with interp2 to look up power for
% measured sea states.
%
% Parameters
% ------------
% device : string
%   Reference model name: "RM3" (wave point absorber), "RM5" (oscillating
%   surge flap), or "RM6" (oscillating water column)
%
% Returns
% ---------
% M : struct
%   M.values : matrix [kW]
%     Device power, Hs bins down the rows and Te bins across the columns
%   M.stat : char
%     'power'
%   M.x_bins : row vector [s]
%     Te bin centers, 0.5 to 20.5 s by 1 s
%   M.y_bins : row vector [m]
%     Hs bin centers, 0.25 to 9.75 m by 0.5 m
%   M.x_edges : row vector [s]
%     Te bin edges
%   M.y_edges : row vector [m]
%     Hs bin edges
%   M.device : char
%     Device name
%   M.description : char
%     Device description from the SAM library
%   M.source : char
%     Data source and license
%
% Examples
% --------
%     >> M = wec_reference_power_matrix("RM3");
%     >> size(M.values)
%     ans =
%         20    21
%
%     >> max(M.values(:))
%     ans =
%        286
%
%     >> % Power [W] for measured sea states, zero outside the matrix
%     >> P = 1000 * interp2(M.x_bins, M.y_bins, M.values, Te, Hm0, 'linear', 0);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    device {mustBeTextScalar, mustBeMember(device, ["RM3", "RM5", "RM6"])}
end

arguments (Output)
    M (1,1) struct
end

file = fullfile(fileparts(mfilename('fullpath')), 'wec_reference_power_matrices', char(device) + ".csv");

% Header comment lines hold the device description and the data source
header = strsplit(strtrim(fileread(file)), newline);
header = header(startsWith(header, '#'));
grid = readmatrix(file, 'CommentStyle', '#');

M = struct();
M.values = grid(2:end, 2:end);
M.stat = 'power';
M.x_bins = grid(1, 2:end);
M.y_bins = grid(2:end, 1)';
M.x_edges = mhkit_define_bins_1d(struct('centers', M.x_bins)).edges';
M.y_edges = mhkit_define_bins_1d(struct('centers', M.y_bins)).edges';
M.device = char(device);
M.description = strtrim(erase(header{1}, '#'));
M.source = strtrim(erase(header{2}, '# Source:'));

end
