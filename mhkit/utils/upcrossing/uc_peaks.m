function peaks = uc_peaks(t, data, inds)
    % Finds the peaks between zero crossings.
    %
    % Parameters:
    % t - Time array.
    % data - Signal time-series.
    % inds - Optional indices for the upcrossing.
    %
    % Returns:
    % peaks - Peak values of the time-series

    if nargin < 3
        inds = upcrossing(t, data);
    end

    peaks = uc_apply_(t, data, @(ind1, ind2) max(data(ind1:ind2)), inds);
end

