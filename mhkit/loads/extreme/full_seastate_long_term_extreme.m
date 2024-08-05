function lte = full_seastate_long_term_extreme(ste_all, weights, x, method, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Return the long-term extreme distribution of a response of
%     interest using the full sea state approach.
%
%     Parameters
%     ----------
%         ste_all : cell array
%             Short-term extreme distribution of the quantity of interest for
%             each sample sea state.
%         weights : array
%             The weights from the full sea state sampling
%
%      Returns
%      -------
%         lte : array or double
%             Long term extreme disribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    ste_all
    weights {mustBeNumeric}
    x
    method
    options.output_py=false;
end

py.importlib.import_module('mhkit');

weights = py.numpy.array(weights);
ste = py.list(ste_all);

result = py.mhkit.loads.extreme.full_seastate_long_term_extreme(ste, weights);

if options.output_py
    lte = result;
else
    stat = py.mhkit_python_utils.scipy_stats.convert_to_array(result, method, x);
    lte = double(stat);
end
