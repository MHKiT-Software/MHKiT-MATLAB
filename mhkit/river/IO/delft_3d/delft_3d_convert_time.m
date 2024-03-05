function result = delft_3d_convert_time(delft_3d_py_object, seconds_run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object.
%    seconds_run: number
%       Time index
%
% Returns
% ---------
%    result: number
%        The closest time index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:get_delft_3d_keys:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');
    end

    python_result = py.mhkit_python_utils.delft_3d_helper.call_convert_time(delft_3d_py_object, seconds_run);

    % Return type can be an numpy int64 or a double. Converting both of these types into MATLAB is
    % more consistent using double conversion
    result = double(python_result) + 1;
end
