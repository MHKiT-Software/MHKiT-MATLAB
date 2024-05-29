function result = delft_3d_index_to_seconds(delft_3d_py_object, index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns the time in seconds corresponding to the given index in the Delft3D dataset.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object, created with `delft_3d_open_netcdf`
%    index: numeric
%       The index for which to retrieve the time in seconds.
%
% Returns
% ---------
%    result: numeric
%        The time in seconds.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_index_to_seconds:InvalidInput', ...
              'Invalid input Delft3D data type: `delft_3d_index_to_seconds` expects a `py.netCDF4._netCDF4.Dataset` object. Please use the `delft_3d_open_netcdf` function to convert Delft3D netCDF files for use with this function.');
    end

    if ~isnumeric(index)
        error('MATLAB:delft_3d_index_to_seconds:InvalidInput', 'Index must be numeric.');
    end

    % Convert to a zero index array for Python
    python_result = py.mhkit.river.io.d3d.index_to_seconds(delft_3d_py_object, int64(index - 1));

    result = double(python_result);
end
