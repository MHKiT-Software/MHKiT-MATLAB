function ste = ste_block_maxima_gev(block_maxima, x, method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Approximate the short-term extreme distribution using the block 
%     maxima method and the Generalized Extreme Value distribution.
%
%     To learn more about the methods that can be used in this function,
%     refer to the scipy.stats.rv_continuous methods documentation!
%
%     Parameters
%     ----------
%         block_maxima : array
%             Block maxima (i.e. largest peak in each block).
%         x : array
%             Input for the statistical function/method
%         method : str
%             Statistical method to apply to resulting data. Options to
%             choose from are: "pdf", "cdf", "ppf", or "sf"
%
%      Returns
%      -------
%         block_maxima : array
%             Block maxima (i.e. largest peak in each block)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(block_maxima,'numeric')
    error('ERROR: block_maxima must be a double array')
end

if ~isa(method,'string')
    error('ERROR: method must be a string')
end

py.importlib.import_module('mhkit');

block_maxima = py.numpy.array(block_maxima);
x = py.numpy.array(x);

result = py.mhkit.loads.extreme.ste_block_maxima_gev(block_maxima);

stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);

ste = double(stat);

end