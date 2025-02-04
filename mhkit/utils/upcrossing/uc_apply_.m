function vals = uc_apply_(t, data, f, inds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply a function over defined intervals in time series data.
%
% Parameters
% ------------
%     t: array
%         Array of time values.
%     data: array
%         Array of data values.
%     f: function handle
%         Function that takes two indices (start, end) and returns a scalar value.
%     inds: array, optional
%         Indices array defining the intervals. If not provided, intervals will be
%         computed using the upcrossing function.
%
% Returns
% ---------
%     vals: array
%         Array of values resulting from applying the function over the defined intervals.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4           % nargin: returns the number of function input arguments given in the call
    % If inds is not provided, compute using upcrossing
    inds = upcrossing(t, data);
end

n = numel(inds) - 1;    % Number of intervals
vals = NaN(1, n);       % Initialize the output array

for i = 1:n
    vals(i) = f(inds(i), inds(i+1));  % Apply the function to each pair of indices
end

end
