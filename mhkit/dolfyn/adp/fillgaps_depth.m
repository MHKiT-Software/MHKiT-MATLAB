function out = fillgaps_depth(var, options)
% Fill gaps (NaN values) in var along depth profile using the specified
% method in options.method
% 
% Parameters
% ----------
% var : struct (from a fieldname in dataset in create_dataset function)
%   The variable to clean
% method : char
%   Interpolation method to use. Default is 'spline'
%
% Returns
% -------
% out : struct (from a fieldname in dataset in create_dataset function)
%  The input struct 'var' with gaps in 'var' interpolated across depth
%
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

    range_dim = 0;
    % find first range dim in var.dims
    for k=1:numel(var.dims) 
       if contains(var.dims{k}, 'range')
           % as soon as range_dim is set once, break out of loop
           range_dim = k;
           break
       end
    end

    if range_dim == 0
        % this means prev loop didnt find 'range' in any dims -> throw error
        error('No range dimension found')
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
    % TODO: interploate over just range dim
    out.data(nans) = interp1(t(not_nans), var.data(not_nans), t(nans), options.method);
end
