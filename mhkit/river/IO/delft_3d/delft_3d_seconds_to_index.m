function result = delft_3d_seconds_to_index(delft_3d_py_object, seconds_elapsed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns the index corresponding to the given number of seconds elapsed in the Delft3D dataset.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object, created with `delft_3d_open_netcdf`
%    seconds_elapsed: numeric
%       The number of seconds elapsed.
%
% Returns
% ---------
%    result: numeric
%        The calculated index.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_seconds_to_index:InvalidInput', ...
              'Invalid input Delft3D data type: `delft_3d_calculate_turbulent_intensity` expects a `py.netCDF4._netCDF4.Dataset` object. Please use the `delft_3d_open_netcdf` function to convert Delft3D netCDF files for use with this function.');
    end

    if ~isnumeric(seconds_elapsed)
        error('MATLAB:delft_3d_seconds_to_index:InvalidInput', 'Seconds elapsed must be numeric.');
    end

    python_result = py.mhkit.river.io.d3d.seconds_to_index(delft_3d_py_object, seconds_elapsed);

    % Account for MATLAB array index difference
    % Note: this changes the index for use in MATLAB, but may cause unexpected behavior in python
    result = double(python_result) + 1;
end
