function n_st = number_of_short_term_peaks(n, t, t_st)
    % Estimate the number of peaks in a specified period.
    % 
    % Parameters
    % ----------
    % n : int
    %     Number of peaks in analyzed timeseries.
    % t : float
    %     Length of time of analyzed timeseries.
    % t_st: float
    %     Short-term period for which to estimate the number of peaks.
    % 
    % Returns
    % -------
    % n_st : float
    %     Number of peaks in short term period.

    if ~isinteger(n) && ~isfloat(n)
        ME = MException('MATLAB:extreme:number_of_short_term_peaks',...
            'n must be an integer');
        throw(ME);
    end
    if ~isinteger(t) && ~isfloat(t)
        ME = MException('MATLAB:extreme:number_of_short_term_peaks',...
            't must be a float');
        throw(ME);
    end
    if ~isinteger(t_st) && ~isfloat(t_st)
        ME = MException('MATLAB:extreme:number_of_short_term_peaks',...
            't_st must be a float');
        throw(ME);
    end

    n_st = n * t_st / t;

end

