function [data, was_row] = mhkit_standardize_user_input_to_column_vectors(data, function_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Standardizes numeric vector input to MHKiT's column-vector convention
%
% MHKiT-MATLAB functions operate on column vectors (Nx1) internally, and
% on matrices (NxK) representing K column-oriented vectors (e.g. K wave
% spectra of length N, one per column). A caller may pass a row vector
% (1xN) instead; this function transposes it to a column and reports the
% original orientation so the caller's output can later be restored via
% mhkit_restore_column_vectors_to_user_input. Matrices (MxN, M>1 and
% N>1) are already column-oriented by convention and pass through
% unchanged.
%
% Parameters
% ------------
% data : numeric vector or matrix
%   User-supplied numeric input to standardize
% function_name : string
%   Name of the calling function, used to scope any error this function
%   raises (e.g. "MHKiT:<function_name>:InvalidInput")
%
% Returns
% ---------
% data : numeric vector or matrix
%   Column-oriented data (Nx1 vector, or an unchanged MxN matrix)
% was_row : logical
%   True if the original input was a row vector (1xN, N>1) and was
%   transposed to a column. Pass this to
%   mhkit_restore_column_vectors_to_user_input to restore the original
%   orientation on shape-preserving outputs only. Reduction outputs
%   (fewer elements than the input, e.g. a moment or a peak period) always
%   stay column-oriented per MHKiT's convention and must not be restored.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    data {mustBeNumeric}
    function_name (1,1) string
end

arguments (Output)
    data
    was_row (1,1) logical
end

if isempty(data)
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        '%s requires non-empty numeric input.', function_name);
end

if ndims(data) > 2
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        ['%s requires 2-D numeric input (a vector or a matrix of column ' ...
         'vectors), got %d-D input.'], function_name, ndims(data));
end

was_row = isrow(data) && ~isscalar(data);

if was_row
    data = data.';
end

end
