function ryv = return_year_value(py_dist, return_year, short_term_period_hr)

%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculate the value from a given distribution corresponding to a particular
%     return year.
% 
%     Parameters
%     ----------
%     py_dist: python callable function object of 1 argument
%         Percentage Point Function (inverse CDF) of short term distribution.
%     return_year: double
%         Return period in years.
%     short_term_period_hr: double
%         Short term period the distribution is created from in hours.
% 
%     Returns
%     -------
%     value: float
%         The value corresponding to the return period from the distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    py_dist
    return_year {mustBeNumeric}
    short_term_period_hr {mustBeNumeric}
end

arguments (Output)
    ryv {mustBeNumeric}
end

% check python_ppf input
mustBeA(py_dist, "py.scipy.stats._distn_infrastructure.rv_continuous_frozen")

% calculate probability of exceedance
pe = 1 - (1 / (return_year * 365.25 * 24 / short_term_period_hr));

ryv = double(py.mhkit_python_utils.scipy_stats.convert_to_array(py_dist, 'ppf', pe));

end