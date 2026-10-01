function F = exceedance_probability(Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the exceedance probability
%
% Tied values are assigned their average rank, matching xarray's rank
% used by MHKiT-Python.
%
% Parameters
% ------------
% Q : struct
%   Discharge data
%     Q.Discharge : vector or matrix [m^3/s]
%       Discharge, one timeseries per column
%     Q.time : vector [datetime or s]
%       Time
%
% Returns
% ---------
% F : struct
%   Exceedance probability data
%     F.F : vector or matrix [%]
%       Exceedance probability, one value per discharge value
%     F.time : vector [s]
%       Time, with datetime converted to epoch seconds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    Q struct
end

arguments (Output)
    F struct
end

time = Q.time;
if any(isdatetime(time))
    time = posixtime(time);
end

% A row vector is a single timeseries, rank it as a column
[discharge, was_row] = mhkit_standardize_user_input_to_column_vectors( ...
    Q.Discharge, 'function_name', mfilename);

n = size(discharge, 1);

rank_ascending = zeros(size(discharge));
for col = 1:size(discharge, 2)
    rank_ascending(:, col) = local_average_rank(discharge(:, col));
end

% Convert to descending rank so the smallest value has the highest
% exceedance probability
rank_descending = n - rank_ascending + 1;
exceedance = 100 * rank_descending / (n + 1);

F.F = mhkit_restore_column_vectors_to_user_input(exceedance, was_row);
F.time = time;

end

function r = local_average_rank(x)
% Assigns ascending ranks (starting at 1) to the elements of x, averaging
% the ranks of tied values. This mirrors the "average" tie-breaking method
% used by xarray/scipy when computing exceedance probability.
x = x(:);
n = numel(x);
[sorted_x, order] = sort(x);
r = zeros(n, 1);

i = 1;
while i <= n
    j = i;
    while j < n && sorted_x(j + 1) == sorted_x(i)
        j = j + 1;
    end
    r(order(i:j)) = (i + j) / 2;
    i = j + 1;
end
end
