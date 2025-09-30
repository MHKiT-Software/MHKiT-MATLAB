function data_out = subtract_mean_from_dimension(data_in, dim)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Subtract mean along specified dimension
%
% Parameters
% ------------
%   data_in : array
%       Input data array
%   dim : double
%       Dimension index along which to subtract mean (typically time dimension)
%
% Returns
% ---------
%   data_out : array
%       Same dimensions as input data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data_in {mustBeNumeric}
    dim (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
end

% Validate dimension is within data array dimensions
if dim > ndims(data_in)
    error('MHKiT:subtract_mean_from_dimension:InvalidDimension', ...
          'Dimension %d exceeds data dimensionality (%d)', dim, ndims(data_in));
end

% Validate data is not empty
if isempty(data_in)
    error('MHKiT:subtract_mean_from_dimension:EmptyData', ...
          'Input data array cannot be empty');
end

data_mean = mean(data_in, dim, 'omitnan');
data_out = data_in - data_mean;
end
