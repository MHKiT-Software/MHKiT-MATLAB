function [pks, result] = peaks_distribution_peaks_over_threshold(peaks, x, method, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Estimate the peaks distribution by fitting a Weibull
%     distribution to the peaks of the response.
%
%     To learn more about the methods that can be used in this function,
%     refer to the scipy.stats.rv_continuous methods documentation!
%
%     Parameters
%     ----------
%         peaks : array
%             Global peaks
%         x : array or double
%             Input for the statistical function/method
%         method : str
%             Statistical method to apply to resulting data. Options to
%             choose from are: "pdf", "cdf", "ppf", or "sf"
%         threshold : double or int (optional)
%             Threshold value. Only peaks above this value will be used.
%             Default value calculated as: `mean(x) + 1.4 * std(x)`
%      Returns
%      -------
%         p : array or double
%             Probability distribution of the peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    peaks
    x
    method
    options.threshold double
end

if ~isa(peaks,'numeric')
    error('ERROR: peaks must be a double array')
end
if ~isa(method,'string')
    error('ERROR: method must be a string')
end

py.importlib.import_module('mhkit');

p = py.numpy.array(peaks);
x = py.numpy.array(x);

if exist("threshold", 'var')==1
    result = py.mhkit.loads.extreme.peaks_distribution_peaks_over_threshold(p, options.threshold);
else
    result = py.mhkit.loads.extreme.peaks_distribution_peaks_over_threshold(p);
end
stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);

pks = double(stat);

end