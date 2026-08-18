function mhkit_verify_column_vector_output(data, function_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Verifies a reduction function's output is a column vector
%
% Runtime check for MHKiT-MATLAB's column-vector output convention.
% Throws if data is not a column vector (or scalar).
%
% Parameters
% ------------
% data : numeric
%   Output to verify
% function_name : string
%   Name of the calling function, used to scope the error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    data {mustBeNumeric}
    function_name (1,1) string
end

if ~iscolumn(data)
    error(sprintf('MHKiT:%s:InvalidOutput', function_name), ...
        '%s must return a column vector, got size [%d %d].', ...
        function_name, size(data,1), size(data,2));
end

end
