function ste = short_term_extreme(t, data, t_st, type, x, method, options)

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
%         t : array
%             Time array
%         data : array
%             Response timeseries
%         t_st : double or int
%             Short time period
%         type : string
%             Method for estimating the short-term extreme distribution
%         x : array or double
%             Input for the statistical function/method
%         method : str
%             Statistical method to apply to resulting data. Options to
%             choose from are: "pdf", "cdf", "ppf", or "sf"
%         output_py : bool (optional)
%             Select if you want to return the native python result for use
%             in long term extreme calculations.
%
%      Returns
%      -------
%         ste : array or double
%             Probability distribution of the peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    t
    data
    t_st
    type
    x
    method
    options.output_py=false;
end


if ~isa(t,'numeric')
    error('ERROR: t must be a double array')
end
if ~isa(data,'numeric')
    error('ERROR: data must be a double array')
end
if ~isa(t_st,'numeric')
    error('ERROR: t_st must be a double')
end
if ~isa(type,'string')
    error('ERROR: type must be a string')
end
if ~isa(x,'numeric')
    error('ERROR: x must be a double array')
end
if ~isa(method,'string')
    error('ERROR: method must be a string')
end

py.importlib.import_module('mhkit');

t = py.numpy.array(t);
data = py.numpy.array(data);

result = py.mhkit.loads.extreme.short_term_extreme(t, data, t_st, type);

if options.output_py
    ste = result;
else
    stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);
    ste = double(stat);
end

end