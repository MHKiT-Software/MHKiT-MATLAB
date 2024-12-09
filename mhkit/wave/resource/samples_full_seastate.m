function [h_sample, t_sample, weight_points] = samples_full_seastate(x1, x2, points_per_interval, return_periods, sea_state_duration, method, bin_size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Sample a sea state between contours of specified return periods.
%    
%    This function is used for the full sea state approach for the
%    extreme load. See Coe et al. 2018 for more details. It was
%    originally part of WDRT.
%    
%    Coe, R. G., Michelen, C., Eckert-Gallup, A., &
%    Sallaberry, C. (2018). Full long-term design response analysis of a
%    wave energy converter. Renewable Energy, 116, 356-366.
%
%     Parameters
%     ----------
%         x1 : array
%             Component 1 data.
%         x2 : array
%             Component 2 data.
%         point_per_interval : int
%             Number of sample points to be calculated per contour interval.
%         return_periods : array
%             Vector of return periods that define the contour intervals in 
%             which samples will be taken. Values must be greater than zero 
%             and must be in increasing order.
%         sea_state_duration : double
%             `x1` and `x2` sample rate (seconds)
%         method : str or list
%             Copula method to apply. Currently only 'PCA' is implemented.
%         bin_size : int
%             Number of data points in each bin (250 is recommended). 
%
%      Returns
%      -------
%         Hs_Samples : array
%             Vector of Hs values for each sample point.
%         Te_Samples : array
%             Vector of Te values for each sample point.
%         weight_points : array
%             Vector of probabilistic weights for each sampling point
%             to be used in risk calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(x1,'numeric')
    error('ERROR: x1 must be a double array')
end
if ~isa(x2,'numeric')
    error('ERROR: x2 must be a double array')
end
if ~isa(points_per_interval,'double')
    error('ERROR: points_per_interval must be an integer')
end
if ~isa(return_periods,'numeric')
    error('ERROR: return_periods must be a double array')
end
if ~isa(sea_state_duration,'numeric')
    error('ERROR: sea_state_duration must be a double')
end
if ~isa(method,'string')
    error('ERROR: method must be a string')
end
if ~isa(bin_size,'double')
    error('ERROR: bin_size must be an integer')
end

x1 = py.numpy.array(x1);
x2 = py.numpy.array(x2);
points_per_interval = int32(points_per_interval);
return_periods = py.numpy.array(return_periods);

result = py.mhkit.wave.contours.samples_full_seastate(x1, x2, points_per_interval, return_periods, sea_state_duration);
result = cell(result);
h_sample = double(result{1});
t_sample = double(result{2});
weight_points = double(result{3});

end

