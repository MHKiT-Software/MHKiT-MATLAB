function result_struct = delft_3d_calculate_turbulent_intensity(delft_3d_py_object, points, intermediate_values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns the turbulent intensity percentage for a given data set for the specified points.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object.
%    points: struct or string or DataFrame
%       Points to interpolate data onto.
%            'cells': interpolates all data onto velocity coordinate system (Default).
%            'faces': interpolates all data onto the TKE coordinate system.
%            DataFrame of x, y, and z coordinates: Interpolates data onto user
%            provided points.
%    intermediate_values: bool (optional)
%       If false, the function will return position and turbulent intensity values.
%       If true, the function will return position (x, y, z) and values needed to calculate
%       turbulent intensity (ucx, ucy, uxz, and turkin1) in a DataFrame. Default is false.
%
% Returns
% ---------
%    result_struct: struct
%        A struct containing the calculated turbulent intensity data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_calculate_turbulent_intensity:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');
    end

    if ~(isstring(points) || isa(points, 'py.pandas.core.frame.DataFrame') || isstruct(points))
        error('MATLAB:delft_3d_calculate_turbulent_intensity:InvalidInput', 'Points must be a string, a struct, or a DataFrame.');
    end

    if ~islogical(intermediate_values)
        error('MATLAB:delft_3d_calculate_turbulent_intensity:InvalidInput', 'Intermediate values must be a boolean.');
    end

    % Convert struct to DataFrame if it's a struct
    if isstruct(points)
        points = convert_numeric_struct_to_dataframe(points);
    end

    % Call Python function to calculate turbulent intensity
    python_result = py.mhkit.river.io.d3d.turbulent_intensity(delft_3d_py_object, pyargs('points', points,'intermediate_values', intermediate_values));

    % Convert Python dataframe to struct
    df = python_result;
    result_struct = convert_numeric_dataframe_to_struct(df);
end
