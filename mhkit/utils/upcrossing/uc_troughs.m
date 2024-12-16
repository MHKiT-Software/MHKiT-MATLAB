function troughs = uc_troughs(t, data, inds)
    % Finds the troughs between zero crossings.
    %
    % Parameters:
    % t - Time array.
    % data - Signal time-series.
    % inds - Optional indices for the upcrossing.
    %
    % Returns:
    % troughs - Trough values of the time-series

    if nargin < 3
        inds = upcrossing(t, data);
    end

    troughs = uc_apply_(t, data, @(ind1, ind2) min(data(ind1:ind2)), inds);
end


