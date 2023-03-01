function out = fillgaps_time(var, options)
% Fill gaps (nan values) in var across time using the specified method in
% options.method
% 
% Parameters
% ----------
% var : struct (from a fieldname in dataset in create_dataset function)
%   The variable to clean
% method : char
%   Interpolation method to use. Default is 'spline'
% Returns
% -------
% out : struct (from a fieldname in dataset in create_dataset function)
%  The input struct 'var' with gaps in 'var' interpolated across time
% See Also
% --------
% Methodology for interpolating across nans from 
% https://www.mathworks.com/matlabcentral/answers/34346-interpolating-nan-s
% create_dataset.m function

% TODO: need to assert that var is data struct from create_dataset here
% set default value for options.method to 'spline'
    arguments
        var;
        options.method char = 'spline';
    end

    time_dim = 0;
    % find first time dim in var.dims
    for k=1:numel(var.dims) 
       if contains(var.dims{k}, 'time')
           % as soon as time_dim is set once, break out of loop
           time_dim = k;
           break
       end
    end

    if time_dim == 0
        % this means prev loop didnt find 'time' in any dims -> throw error
        error('No time dimension found')
    end

    % now interpolate the nans away

    % get array of true false for which are nans in var.data
    % not_nans has 0 if nan and 1 if number (since ~)
    nans = isnan(var.data);
    not_nans = ~nans;
    % t is the array {1,2,3,....,num_elements_in_var.data}
    t = 1:numel(var.data); 
    % output is var except the nans have been interpolated:
    out = var;
    % reassign any nans in out.data to be interpolated values at the nan
    % positions
    out.data(nans) = interp1(t(not_nans), var.data(not_nans), t(nans), options.method);

