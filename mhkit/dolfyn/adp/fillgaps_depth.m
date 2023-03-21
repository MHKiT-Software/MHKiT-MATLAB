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
%   see fillmissing method for other options 
%   https://www.mathworks.com/help/matlab/ref/fillmissing.html#bva1z1c-method
% options.maxgap : double 
%   Maximum gap of missing data to interpolate across. Default is None
%
% Returns
% -------
% out : struct (from a fieldname in dataset in create_dataset function)
%  The input struct 'var' with gaps in 'var' interpolated across depth
%
% See Also
% --------
% create_dataset.m function

% TODO: need to assert that var is data struct from create_dataset here
% set default value for options.method to 'spline'
    arguments
        var;
        options.method char = 'spline';
        options.maxgap double = NaN;
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

    % now interpolate the nans away, note fillmissing will look for NaNs in
    % var.data as the missing entries
    out = var;
    if isnan(options.maxgap) %if maxgap is Nan, no maxgap was set:
        out.data = fillmissing(var.data,options.method,range_dim);
    end
    if ~isnan(options.maxgap) % if maxgap is not NaN, a maxgap is set:
        out.data = fillmissing(var.data,options.method,range_dim,MaxGap=options.maxgap);
    end
    end
