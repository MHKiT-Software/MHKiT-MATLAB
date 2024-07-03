function ste = ste_peaks(peaks_distribution, npeaks, x, method)

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
%         peaks_distribution : scipy.stats.rv_frozen
%             Probability distribution of the peaks
%         npeaks : double or int
%             Number of peaks in short term period
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

if ~isa(npeaks,'numeric')
    error('ERROR: npeaks must be a double')
end
if ~isa(method,'string')
    error('ERROR: method must be a string')
end

py.importlib.import_module('mhkit');

x = py.numpy.array(x);

result = py.mhkit.loads.extreme.ste_peaks(peaks_distribution, npeaks);

stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);

ste = double(stat);

end