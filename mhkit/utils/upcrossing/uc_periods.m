function periods = uc_periods(t, data, inds)
    % Calculates the period between zero crossings.
    %
    % Parameters:
    % t - Time array.
    % data - Signal time-series.
    % inds - Optional indices for the upcrossing.
    %
    % Returns:
    % periods - Period values of the time-series

    if nargin < 3
        inds = upcrossing(t, data);
    end

    periods = uc_apply_(t, data, @(ind1, ind2) t(ind2) - t(ind1), inds);
end