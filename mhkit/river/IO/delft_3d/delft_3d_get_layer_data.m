function result_struct = delft_3d_get_layer_data(delft_3d_py_object, variable, layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns layer data from a Delft 3D netCDF object.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object.
%    variable: string
%       Variable in the Delft3D output
%    layer: int
%       Layer index
%
% Returns
% ---------
%    result_struct: struct
%        A struct containing the layer data with fields corresponding to the
%        DataFrame columns, including the DataFrame itself stored as 'df'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Guard to check type of delft_3d_py_object
    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_get_layer_data:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');
    end

    % Guard to check type of variable
    if ~isstring(variable)
        error('MATLAB:delft_3d_get_layer_data:InvalidInput', 'Variable must be a string.');
    end

    % Guard to check type of layer
    if ~isnumeric(layer) || ~isscalar(layer) || layer ~= fix(layer)
        error('MATLAB:delft_3d_get_layer_data:InvalidInput', 'Layer must be a scalar integer.');
    end


    % Call Python function to get layer data
    python_result = py.mhkit.river.io.d3d.get_layer_data(delft_3d_py_object, variable, py.int(int64(layer)));

    % Convert Python DataFrame to struct with vectorized conversion
    result_struct = convert_numeric_dataframe_to_struct(python_result);
end
