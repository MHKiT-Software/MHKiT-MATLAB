function custom_vals = uc_custom(t, data, func, inds)
% Applies a custom function to the time-series data between upcrossing points.
%
% Parameters:
%------------
%     t: array
%           Array of time values.
%     data: array
%         Array of data values.
%     func: function handle
%         Function to apply between the zero crossing periods.
%     inds: Array, optional
%         Indices for the upcrossing.
%
% Returns:
% ---------
%     custom_vals: array
%           Custom values of the time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    inds = upcrossing(t, data);
end
if ~isa(func, 'function_handle')
    error('func must be a function handle');
end

custom_vals = uc_apply_(t, data, func, inds);
end
