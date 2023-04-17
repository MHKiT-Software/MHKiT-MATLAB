function out = fillgaps_time(var, options)
% Fill gaps (NaN values) in var across time using the specified method in
% options.method
% 
% Parameters
% ----------
% var : struct (from a fieldname in dataset in create_dataset function)
%   The variable to clean
% options.method : char
%   Interpolation method to use. Default is 'spline'
%   see fillmissing method for other options 
%   https://www.mathworks.com/help/matlab/ref/fillmissing.html#bva1z1c-method
% options.maxgap : double 
%   Maximum gap of missing data to interpolate across. Default is None
%
% Returns
% -------
% out : struct (from a fieldname in dataset in create_dataset function)
%  The input struct 'var' with gaps in 'var' interpolated across time
%
% See Also
% --------
% create_dataset.m function

% set default value for options.method to 'spline'
    arguments
        var;
        options.method char = 'spline';
        options.maxgap double = NaN;
    end
    % Make sure that var is a structure
    if ~isstruct(var)
        ME = MException('MATLAB:fillgaps_time','var must be a structure');
        throw(ME);
    end

    % Make sure that ds contains the dolyn fields
    if ~isfield(var,'coords') || ~isfield(var,'data') 
        ME = MException('MATLAB:fillgaps_time',['The provided data ' ...
            'structure does not appear to have been created by dolfyn']);
        throw(ME);
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

    % if time_dim is still 0, then no time_dim was found in var.dims
    if time_dim == 0
        ME = MException('MATLAB:fillgaps_time',['No time dimension found']);
        throw(ME);
    end

    % now interpolate the nans away
    out = var;
    if isnan(options.maxgap) %if maxgap is Nan, no maxgap was set:
        out.data = fillmissing(var.data,options.method,time_dim);
    end
    if ~isnan(options.maxgap) % if maxgap is not NaN, a maxgap is set:
        out.data = fillmissing(var.data,options.method,time_dim,'MaxGap', options.maxgap);
    end
end
