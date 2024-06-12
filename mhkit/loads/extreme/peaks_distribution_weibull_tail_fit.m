function [pks, result] = peaks_distribution_weibull_tail_fit(peaks, x, method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Estimate the peaks distribution by using the Weibull tail fit method.
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
%
%      Returns
%      -------
%         p : array or double
%             Probability distribution of the peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(peaks,'numeric')
    error('ERROR: peaks must be a double array')
end

if ~isa(method,'string')
    error('ERROR: method must be a string')
end

py.importlib.import_module('mhkit');

p = py.numpy.array(peaks);
x = py.numpy.array(x);

result = py.mhkit.loads.extreme.peaks_distribution_weibull_tail_fit(p);

stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);

pks = double(stat);

end