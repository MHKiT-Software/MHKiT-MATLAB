function mhkit_verify_is_column_vector(data, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Verifies a value is a column vector (or scalar), input or output
%
% Runtime sanity check for MHKiT-MATLAB's column-vector convention. Use
% on a spectral statistic function's output, or on an intermediate value
% (e.g. frequency, time) to catch a caller's mistake early rather than
% let it fail confusingly deeper in the computation. Throws if data is
% not a column vector or scalar.
%
% Parameters
% ------------
% data : numeric, datetime, or duration
%   Value to verify
% function_name : string (optional)
%   Name of the calling function, used to scope the error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    data
    options.function_name (1,1) string = mfilename
end

if ~(isnumeric(data) || isdatetime(data) || isduration(data))
    error(sprintf('MHKiT:%s:InvalidInput', options.function_name), ...
        '%s expected numeric, datetime, or duration, got %s.', options.function_name, class(data));
end

if ~iscolumn(data)
    error(sprintf('MHKiT:%s:InvalidOutput', options.function_name), ...
        '%s expected a column vector, got size [%d %d].', ...
        options.function_name, size(data,1), size(data,2));
end

end
