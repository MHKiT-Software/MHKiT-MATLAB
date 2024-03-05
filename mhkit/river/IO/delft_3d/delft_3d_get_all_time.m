function result = delft_3d_get_all_time(delft_3d_py_object)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns all time values from a Delft 3D netCDF object.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object.
%
% Returns
% ---------
%    result: double
%        An array containing all time values from the Delft 3D netCDF object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_get_all_time:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');
    end

    python_result = py.mhkit.river.io.d3d.get_all_time(delft_3d_py_object);

    result = double(python_result);
end
