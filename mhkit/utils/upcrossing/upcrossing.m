function zeroUpCrossings_index = upcrossing(t, data)
    % Finds the zero upcrossing points.
    %
    % Parameters:
    % t - Time array.
    % data - Signal time series.
    %
    % Returns:
    % Zero crossing indices

    if ~isvector(t) || ~isvector(data)
        error('t and data must be 1D arrays');
    end
    
    % eliminate zeros
    zeroMask = (data == 0);
    data(zeroMask) = 0.5 * min(abs(data));

    % zero up-crossings
    diff_data = diff(sign(data));
    zeroUpCrossings_mask = (diff_data == 2) | (diff_data == 1);
    zeroUpCrossings_index = find(zeroUpCrossings_mask);
end

