function result_struct = delft_3d_calculate_variable_interpolation(delft_3d_py_object, variables, points, edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interpolates multiple variables from the Delft3D onto the same points.
%
% Parameters
% ------------
%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset
%       A netCDF python object.
%    variables: cell array of strings
%        Name of variables to interpolate, e.g. 'turkin1', 'ucx', 'ucy', and 'ucz'.
%        The full list can be found using "data.variables.keys()" in the console.
%    points: string or struct or DataFrame
%        The points to interpolate data onto.
%          'cells': interpolates all data onto the Delft3D cell coordinate system (Default)
%          'faces': interpolates all data onto the Delft3D face coordinate system
%          DataFrame of x, y, and waterdepth coordinates: Interpolates data onto user
%          provided points. Can be created with `create_points` function.
%    edges: string, optional
%        If edges is set to 'nearest', the code will fill in NaN values with nearest
%        interpolation. Otherwise, only linear interpolation will be used.
%
% Returns
% ---------
%    result_struct: struct
%        A struct containing the interpolated variables on specified grid points saved
%        under the input variable names and the x, y, and waterdepth coordinates of those points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')
        error('MATLAB:delft_3d_calculate_variable_interpolation:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');
    end

    if ~iscell(variables) || ~all(cellfun(@ischar, variables))
        error('MATLAB:delft_3d_calculate_variable_interpolation:InvalidInput', 'Variables must be a cell array of strings.');
    end

    if ~(isstring(points) || isa(points, 'py.pandas.core.frame.DataFrame') || isstruct(points))
        error('MATLAB:delft_3d_calculate_variable_interpolation:InvalidInput', 'Points must be a string or a DataFrame.');
    end

    if nargin < 4
        edges = '';
    end

    if ~isstring(edges) && ~isempty(edges)
        error('MATLAB:delft_3d_calculate_variable_interpolation:InvalidInput', 'Edges must be a string or empty.');
    end

    % Convert struct to DataFrame if it's a struct
    if isstruct(points)
        points = convert_numeric_struct_to_dataframe(points);
    end

    % Call Python function to perform variable interpolation
    python_result = py.mhkit.river.io.d3d.variable_interpolation(delft_3d_py_object, py.list(variables), pyargs('points', points, 'edges', edges));

    % Convert Python dataframe to struct
    df = python_result;
    result_struct = convert_numeric_dataframe_to_struct(df);
end
