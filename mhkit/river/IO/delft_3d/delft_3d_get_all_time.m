function result = delft_3d_get_all_time(delft_3d_py_object)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns all time values from a Delft 3D netCDF object.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object, created with `delft_3d_open_netcdf`
%
% Returns
% ---------
%    result: double
%        An array containing all time values from the Delft 3D netCDF object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
    error('MATLAB:delft_3d_get_all_time:InvalidInput', ...
        'Invalid input Delft3D data type: `delft_3d_get_all_time` expects a `py.netCDF4._netCDF4.Dataset` object. Please use the `delft_3d_open_netcdf` function to convert Delft3D netCDF files for use with this function.');
end

    python_result = py.mhkit.river.io.d3d.get_all_time(delft_3d_py_object);

    % Handle masked arrays using numpy.ma.getdata
    if isa(python_result, 'py.numpy.ma.MaskedArray')
        % Extract the underlying data
        data_array = py.numpy.ma.getdata(python_result);
    else
        % Convert directly if not masked
        data_array = py.numpy.array(python_result);
    end

    % Ensure data is of type float
    float_array = data_array.astype('float');

    result = double(float_array);
end
