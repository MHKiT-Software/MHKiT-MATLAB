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

    python_result = py.list(python_result);

disp(python_result);
disp(class(python_result));

result = double(python_result);
end
