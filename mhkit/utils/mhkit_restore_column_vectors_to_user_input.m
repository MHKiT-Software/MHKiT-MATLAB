function data = mhkit_restore_column_vectors_to_user_input(data, was_row)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Restores column-oriented data to the caller's original vector orientation
%
% Sister function to mhkit_standardize_user_input_to_column_vectors. Only
% use this for shape-preserving functions, where the output has the same
% number of elements as the standardized input (e.g. surface_elevation).
% Reduction functions (fewer output elements than input, e.g.
% frequency_moment, peak_period) always emit column-oriented output per
% MHKiT's vector convention and must not call this function.
%
% Parameters
% ------------
% data : numeric vector or matrix
%   Column-oriented output to restore
% was_row : logical
%   Orientation flag returned by
%   mhkit_standardize_user_input_to_column_vectors
%
% Returns
% ---------
% data : numeric vector or matrix
%   Transposed back to a row vector if was_row is true, unchanged
%   otherwise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    data {mustBeNumeric}
    was_row (1,1) logical
end

arguments (Output)
    data
end

if was_row
    data = data.';
end

end
