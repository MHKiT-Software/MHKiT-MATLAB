function m = mhkit_nanmean(x, dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Compute mean ignoring NaN values
%
%     Computes the mean of input array, ignoring NaN values. This function 
%     provides equivalent functionality to Python numpy.nanmean()%
%
% Parameters
% ------------
%   x : numeric array
%       Input array containing numeric data with potential NaN values
%   dim : positive integer, optional
%       Dimension along which to compute mean. If not specified, 
%       operates along first non-singleton dimension [dimensionless]
%
% Returns
% ---------
%   m : numeric array
%       Mean values with NaN values excluded. Output dimensions depend
%       on input array and specified dimension [same units as input]
%
% Examples
% ---------
%   % Vector example
%   x = [1, 2, NaN, 4, 5];
%   m = mhkit_nanmean(x);  % Returns 3 (mean of [1,2,4,5])
%
%   % Matrix example - column-wise mean
%   x = [1 2 3; NaN 5 6; 7 8 NaN];
%   m = mhkit_nanmean(x, 1);  % Returns [4, 5, 4.5]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    x {mustBeNumeric}
    dim {mustBeNumeric, mustBeInteger} = []
end

% Validate input array is not empty
if isempty(x)
    error('MHKiT:mhkit_nanmean:EmptyInput', ...
        'Input array cannot be empty');
end

% Handle optional dimension parameter
if isempty(dim)
    % Find first non-singleton dimension
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
else
    % Validate provided dimension is positive
    if dim < 1
        error('MHKiT:mhkit_nanmean:InvalidDimension', ...
            'Dimension must be a positive integer, got %d', dim);
    end
end

% Validate dimension is within valid range
if dim > ndims(x)
    error('MHKiT:mhkit_nanmean:InvalidDimension', ...
        'Dimension %d exceeds the number of dimensions in input array (%d)', ...
        dim, ndims(x));
end

% Find non-NaN elements
valid = ~isnan(x);

if ~any(valid(:))
    % All values are NaN
    m = NaN(size(x));
    return;
end

% Sum valid values and count them
x_sum = sum(x .* valid, dim, 'omitnan');
x_count = sum(valid, dim);

% Compute mean, handling case where all values in a slice are NaN
m = x_sum ./ x_count;

% Set result to NaN where no valid values exist
m(x_count == 0) = NaN;

end
