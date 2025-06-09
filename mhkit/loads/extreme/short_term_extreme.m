function ste = short_term_extreme(t, data, t_st, type, x, method, options)
% SHORT_TERM_EXTREME Estimate short-term extreme value statistics
%
%   ste = short_term_extreme(t, data, t_st, type, x, method)
%
%   Estimates extreme value distributions from a timeseries of the response
%   using either block maxima or peak-over-threshold methods.
%
% Parameters
% ----------
%   t       : time vector
%   data    : response signal
%   t_st    : short-term period (scalar)
%   type    : method string:
%             'peaks_weibull', 'peaks_weibull_tail_fit',
%             'peaks_over_threshold', 'block_maxima_gev', 'block_maxima_gumbel'
%   x       : vector of values to evaluate the distribution at
%   method  : 'pdf', 'cdf', 'ppf' (inverse CDF), or 'sf' (1 - CDF)
%   options : struct with optional fields:
%             .threshold (required for peaks_over_threshold)
%
% Returns
% -------
%   ste     : values of the estimated distribution at x

arguments
    t (:,1) double
    data (:,1) double
    t_st (1,1) double
    type (1,:) char
    x (:,1) double
    method (1,:) char
    options.threshold double = NaN
end

% Error checking
if length(t) ~= length(data)
    error('Time and data must be the same length');
end

% Sampling interval and frequency
dt = mean(diff(t));
fs = 1/dt;

switch lower(type)
    case 'peaks_weibull'
        peaks = find_peaks(data);
        nst = round(length(peaks) * t_st / duration);
        pd_st = fit_weibull(peaks, nst, x);

    case 'peaks_weibull_tail_fit'
        peaks = find_peaks(data);
        pd_st = fit_weibull_tail(peaks, x);

    case 'peaks_over_threshold'
        if isnan(options.threshold)
            error('You must provide options.threshold for this method.');
        end
        peaks = data(data > options.threshold);
        pd_st = fit_pareto(peaks, x);

    case 'block_maxima_gev'
        block_size = round(t_st * fs);
        blocks = reshape_truncated(data, block_size);
        maxima = max(blocks, [], 1);
        pd_st = fit_gumbel(maxima, x);

    case 'block_maxima_gumbel'
        block_size = round(t_st * fs);
        blocks = reshape_truncated(data, block_size);
        maxima = max(blocks, [], 1);
        pd_st = fit_gev(maxima, x);

    otherwise
        error('Unknown method: %s', type);
end

switch lower(method)
    case 'pdf'
        ste = pd_st.pdf;
    case 'cdf'
        ste = pd_st.cdf;
    case {'ppf', 'inv'}
        ste = pd_st.ppf;
    case 'sf'
        ste = 1 - pd_st.cdf;
    case 'expect'
        % Numerical expectation: integral of x * pdf(x)
        dx = mean(diff(x));
        ste = sum(x .* pd_st.pdf) * dx;
    otherwise
        error('Invalid method: %s. Choose from "pdf", "cdf", "ppf", "sf", "expect".', method);
end
end


function peaks = find_peaks(data)
    locs = find(data(2:end-1) > data(1:end-2) & data(2:end-1) > data(3:end)) + 1;
    peaks = data(locs);
end

function blocks = reshape_truncated(data, block_size)
    N = floor(length(data)/block_size) * block_size;
    blocks = reshape(data(1:N), block_size, []);
end

function pd = fit_weibull_tail(data, x)
    % Manual peak filtering
    locs = find(data(2:end-1) > data(1:end-2) & data(2:end-1) > data(3:end)) + 1;
    peaks = data(locs);

    % Tail selection
    threshold = prctile(peaks, 95);
    tail_peaks = peaks(peaks > threshold);

    % Grid search over Weibull parameters
    scale_vals = linspace(mean(tail_peaks)*0.5, mean(tail_peaks)*2, 50);
    shape_vals = linspace(0.5, 5, 50);

    min_nll = inf;

    for a = scale_vals
        for b = shape_vals
            pdf_vals = weibull_pdf(tail_peaks, a, b);
            if all(pdf_vals > 0 & isfinite(pdf_vals))
                nll = -sum(log(pdf_vals));
                if nll < min_nll
                    min_nll = nll;
                    best_scale = a;
                    best_shape = b;
                end
            end
        end
    end

    pd.type = 'Weibull';
    pd.scale = best_scale;
    pd.shape = best_shape;
    pd.pdf = weibull_pdf(x, best_scale, best_shape);
    pd.cdf = weibull_cdf(x, best_scale, best_shape);
    pd.ppf = weibull_ppf(x, best_scale, best_shape);
end

function pd = fit_weibull(data, nst, x)
    % Same fitting approach as the tail, just on all peaks
    scale_vals = linspace(mean(data)*0.5, mean(data)*2, 50);
    shape_vals = linspace(0.5, 5, 50);

    min_nll = inf;

    for a = scale_vals
        for b = shape_vals
            pdf_vals = weibull_pdf(data, a, b);
            if all(pdf_vals > 0 & isfinite(pdf_vals))
                nll = -sum(log(pdf_vals));
                if nll < min_nll
                    min_nll = nll;
                    best_scale = a;
                    best_shape = b;
                end
            end
        end
    end

    pd.type = 'Weibull';
    pd.scale = best_scale;
    pd.shape = best_shape;
    pd.pdf = weibull_pdf(x, best_scale, best_shape);
    pd.cdf = weibull_cdf(x, best_scale, best_shape).^nst;  % scaled for ST realization
    pd.ppf = weibull_ppf(x, best_scale, best_shape);
end

function p = weibull_pdf(x, scale, shape)
    p = (shape./scale) .* (x./scale).^(shape - 1) .* exp(-(x./scale).^shape);
end

function c = weibull_cdf(x, scale, shape)
    c = 1 - exp(-(x./scale).^shape);
end

function x_ppf = weibull_ppf(p, scale, shape)
    x_ppf = scale .* (-log(1 - p)).^(1/shape);
end
