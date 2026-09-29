function F=exceedance_probability(Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the exceedance probability
%
% Parameters
% ----------
%     Q : Discharge data [m3/s]
%
%         structure of form:
%
%           Q.Discharge
%
%           Q.time
%
% Returns
% -------
%     F : Structure
%
%
%         F.F: Exceedance probability [unitless]
%
%         F.time: time [epoch time (s)]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

discharge = Q.Discharge;
n = size(discharge, 1);

rank_ascending = zeros(size(discharge));
for col = 1:size(discharge, 2)
    rank_ascending(:, col) = local_average_rank(discharge(:, col));
end

% Convert to descending rank so the smallest value has the highest
% exceedance probability
rank_descending = n - rank_ascending + 1;
F.F = 100 * rank_descending / (n + 1);
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

