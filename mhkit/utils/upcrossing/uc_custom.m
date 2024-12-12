function custom_vals = uc_custom(t, data, func, inds)
    % Applies a custom function to the time-series data between upcrossing points.
    %
    % Parameters:
    % t - Time array.
    % data - Signal time-series.
    % func - Function to apply between the zero crossing periods.
    % inds - Optional indices for the upcrossing.
    %
    % Returns:
    % custom_vals - Custom values of the time-series

    if nargin < 4
        inds = upcrossing(t, data);
    end
    if ~isa(func, 'function_handle')
        error('func must be a function handle');
    end

    custom_vals = uc_apply_(t, data, func, inds);
end
