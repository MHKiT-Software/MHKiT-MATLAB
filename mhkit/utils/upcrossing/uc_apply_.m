function vals = uc_apply_(t, data, f, inds)
    % Apply a function `f` over intervals defined by `inds`. If `inds` is
    % not provided, compute the indices using the upcrossing function.
    
    % Input:
    %   t : Array of time values.
    %   data : Array of data values.
    %   f : Function handle that takes two indices (start, end) and returns a scalar value.
    %   inds : Indices array (optional). If not provided, `upcrossing` will be used.
    
    % Output:
    %   vals : Array of values resulting from applying `f` over the intervals defined by `inds`.

    if nargin < 4           % nargin: returns the number of function input arguments given in the call
        % If inds is not provided, compute using upcrossing
        inds = upcrossing(t, data);
    end
    
    n = numel(inds) - 1;    % Number of intervals
    vals = zeros(1, n);       % Initialize the output array
    
    for i = 1:n
        vals(i) = f(inds(i), inds(i+1));  % Apply the function to each pair of indices
    end
end
