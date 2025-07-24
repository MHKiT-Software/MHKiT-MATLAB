function [best_threshold,best_threshold_unit] = automatic_hs_threshold(peaks, sampling_rate, initial_threshold_range, max_refinement)
%%%%%%%%%%%%%%%%%%%%%
%{
    Find the best significant wave height threshold for the
    peaks-over-threshold method.

    This method was developed by:

    > Neary, V. S., S. Ahn, B. E. Seng, M. N. Allahdadi, T. Wang, Z. Yang and R. He (2020).
    > "Characterization of Extreme Wave Conditions for Wave Energy Converter Design and
    >   Project Risk Assessment.”
    > J. Mar. Sci. Eng. 2020, 8(4), 289; https://doi.org/10.3390/jmse8040289.

    Please cite this paper if using this method.

    After all thresholds in the initial range are evaluated, the search
    range is refined around the optimal point until either (i) there
    is minimal change from the previous refinement results, (ii) the
    number of data points become smaller than about 1 per year, or (iii)
    the maximum number of iterations is reached.

    Parameters
    ----------
    peaks: array
        Peak values of the response time-series.
    sampling_rate: double
        Sampling rate in hours.
    initial_threshold_range: array(double, double, double)
        Initial range of thresholds to search. Described as
        (min, max, step).
    max_refinement: int32
        Maximum number of times to refine the search range.

    Returns
    -------
    array[float, float]
        The best threshold and its corresponding unit.
%}

arguments
    peaks {mustBeNumeric}
    sampling_rate {mustBeNumeric}
    initial_threshold_range {mustBeNumeric} = [0.99, 0.995, 0.001]
    max_refinement {mustBeInteger} = 5
end

py.importlib.import_module('mhkit')

peaks = py.numpy.array(peaks);
initial_threshold_range = py.tuple(initial_threshold_range);

result = py.mhkit.loads.extreme.peaks.automatic_hs_threshold(...
    peaks, sampling_rate, initial_threshold_range, max_refinement);

best_threshold = double(result{1});
best_threshold_unit = double(result{2});

end

