function ste = short_term_extreme(t, data, t_st, method)
    % Approximate the short-term  extreme distribution from a
    % timeseries of the response using chosen method.
    % 
    % The availabe methods are: peaks_weibull, peaks_weibull_tail_fit,
    % peaks_over_threshold, block_maxima_gev, and block_maxima_gumbel.
    % For the block maxima methods the timeseries needs to be many times
    % longer than the short-term period. For the peak-fitting methods the
    % timeseries can be of arbitrary length.
    % 
    % Parameters
    % ----------
    % t: array
    %     Time array.
    % data: array
    %     Response timeseries.
    % t_st: float
    %     Short-term period.
    % method : string
    %     Method for estimating the short-term extreme distribution.
    % 
    % Returns
    % -------
    % ste: peaks_distribution or ste_block_maxima class object
    %         Short-term extreme distribution.

    assert(isvector(t),'t must be array')
    assert(isvector(data),'data must be array')
    assert(isfloat(t_st),'t_st must be of type float')
    assert(isfloat(t_st),'t_st must be of type float')
    if ~ischar(method) && ~isstring(method)
        ME = MException('MATLAB:extreme:short_term_extreme', ...
            "Method must be of type char or string");
        throw(ME);
    end

    peaks_dist = peaks_distribution();
    fit_maxima = ste_block_maxima();

    peaks_methods = struct(  ...
        'peaks_weibull',{@peaks_dist.weibull},...
        'peaks_weibull_tail_fit',{@peaks_dist.weibull_tail_fit},...
        'peaks_over_threshold',{@peaks_dist.peaks_over_threshold});

    blockmaxima_methods = struct(  ...
        'block_maxima_gev',{@fit_maxima.gev},...
        'block_maxima_gumbel',{@fit_maxima.gumbel});

    if isfield(peaks_methods, method)
        [~, peaks] = global_peaks(t, data);
        npeaks = length(peaks);
        time = t(end) - t(1);
        nst = number_of_short_term_peaks(npeaks, time, t_st);        
        feval(peaks_methods.(method),peaks);
        ste = ste_peaks(peaks_dist, nst);
    elseif isfield(blockmaxima_methods, method)
        maxima = block_maxima(t, data, t_st);
        feval(blockmaxima_methods.(method),maxima);
        ste = fit_maxima;
    else
        ME = MException('MATLAB:extreme:short_term_extreme', ...
            "The supplied method: %s is not available. Available methods " + ...
            "are:\n\t-peaks_weibull\n\t-peaks_weibull_tail_fit\n\t" + ...
            "-peaks_over_threshold\n\t-block_maxima_gev\n\t-" + ...
            "block_maxima_gumbel\n", method);
        throw(ME);
    end
    
end

