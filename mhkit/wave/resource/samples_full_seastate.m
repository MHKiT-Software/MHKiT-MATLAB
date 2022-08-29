function [h_sample, t_sample, weight_points] = samples_full_seastate( ...
    x1, x2, points_per_interval, return_periods, sea_state_duration, ...
    options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample a sea state between contours of specified return periods.
% 
% This function is used for the full sea state approach for the
% extreme load. See Coe et al. 2018 for more details. It was
% originally part of WDRT.
% 
% Coe, R. G., Michelen, C., Eckert-Gallup, A., &
% Sallaberry, C. (2018). Full long-term design response analysis of a
% wave energy converter. Renewable Energy, 116, 356-366.
% 
% Parameters
% ----------
% x1: array
%     Component 1 data
% x2: array
%     Component 2 data
% points_per_interval : int
%     Number of sample points to be calculated per contour interval.
% return_periods: array
%     Vector of return periods that define the contour intervals in
%     which samples will be taken. Values must be greater than zero
%     and must be in increasing order.
% sea_state_duration : int or float
%     x1 and x2 sample rate (seconds)
% method: string or list
%     Copula method to apply. Currently only 'PCA' is implemented.
% bin_size : int
%     Number of data points in each bin
% 
% Returns
% -------
% Hs_Samples: array
%     Vector of Hs values for each sample point.
% Te_Samples: array
%     Vector of Te values for each sample point.
% weight_points: array
%     Vector of probabilistic weights for each sampling point
%     to be used in risk calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    x1
    x2
    points_per_interval
    return_periods
    sea_state_duration
    options.method = "PCA";
    options.bin_size = 250;
end

method = options.method;
bin_size = options.bin_size;

if ~strcmpi(method,"PCA")
    ME = MException('MATLAB:wave.resource:samples_full_seastate',...
        ['Full sea state sampling is currently only implemented using ' ...
        'the PCA method.']);
    throw(ME);
end
assert(isvector(x1),'x1 must be array')
assert(isvector(x2),'x2 must be array')
assert(isfloat(points_per_interval) || isinteger(points_per_interval), ...
    'points_per_interval must be of type int')
assert(isvector(return_periods),'return_periods must be array')
assert(isfloat(sea_state_duration) || isinteger(sea_state_duration), ...
    'sea_state_duration must be of type int or float')
assert(ischar(method) || isstring(method), 'method must be of type string')
assert(isfloat(bin_size) || isinteger(bin_size), ...
    'bin_size must be of type int')

bin_size = int32(bin_size);

pca_fit = principal_component_analysis(x1, x2, 'bin_size',bin_size);

h_sample = 0;
t_sample = 0;
weight_points = 0;
end

