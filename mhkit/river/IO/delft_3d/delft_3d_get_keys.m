function result = delft_3d_get_keys(delft_3d_py_object)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns a struct of the key/values of a Delft 3D netCDF object.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object, created with `delft_3d_open_netcdf`
%
% Returns
% ---------
%    result: struct
%        A struct containing the keys and their corresponding values from
%        the Delft 3D netCDF object.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:get_delft_3d_keys:InvalidInput', ...
              'Invalid input Delft3D data type: `delft_3d_get_keys` expects a `py.netCDF4._netCDF4.Dataset` object. Please use the `delft_3d_open_netcdf` function to convert Delft3D netCDF files for use with this function.');
    end

    % get_d3d_keys returns a list with the following elements:
    %   1: List of keys
    %   2: List of values
    python_result = py.mhkit_python_utils.delft_3d_helper.get_d3d_keys(delft_3d_py_object);

    keys = cell(python_result{1});  % Convert python string array to matlab char array
    values = cell(python_result{2});  % Convert python string array to matlab char array

    result = struct();

    for i = 1:numel(keys)
        key = char(keys{i});
        value = char(values{i});
        result.(key) = value;
    end
end
