%function result = delft_3d_interpolate_to_centerline(points, values, xi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Returns
%%
%% Parameters
%% ------------
%%    points: 2-D ndarray of floats
%%       Data point coordinates
%%    values: ndarray of float
%%       Data values
%%    xi: 2-D ndarray of floats
%%       points at which to interoplate data
%%
%% Returns
%% ---------
%%    result: ndarray
%%        Array of interpolated values
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    % Transposed to get the same shape as an input DataFrame
%    points = py.numpy.transpose(py.numpy.asarray(points));
%    values = py.numpy.asarray(values);
%    % Transposed to get the same shape as an input DataFrame
%    xi = py.numpy.transpose(py.numpy.asarray(xi));

%    python_result = py.scipy.interpolate.griddata(points, values, xi);

%    result = double(python_result);
%end

function result = delft_3d_interpolate_to_centerline(points, values, xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolates data to points along a centerline.
%
% Parameters
% ------------
%    points: 2-D ndarray of floats
%       Data point coordinates
%    values: ndarray of float
%       Data values
%    xi: 2-D ndarray of floats
%       Points at which to interpolate data
%
% Returns
% ---------
%    result: ndarray
%        Array of interpolated values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Guard to check type of points
    if ~isnumeric(points)
        error('MATLAB:delft_3d_interpolate_to_centerline:InvalidInput', 'points must be a numeric array.');
    end

    % Guard to check type of values
    if ~isnumeric(values)
        error('MATLAB:delft_3d_interpolate_to_centerline:InvalidInput', 'values must be a numeric array.');
    end

    % Guard to check type of xi
    if ~isnumeric(xi)
        error('MATLAB:delft_3d_interpolate_to_centerline:InvalidInput', 'xi must be a numeric array.');
    end

    % Transposed to get the same shape as an input DataFrame
    % Note that points and xi can be 2d arrays and values can only be a 1d array
    points = py.numpy.transpose(py.numpy.asarray(points));
    values = py.numpy.asarray(values);
    xi = py.numpy.transpose(py.numpy.asarray(xi));

    python_result = py.scipy.interpolate.griddata(points, values, xi);

    result = double(python_result);
end
