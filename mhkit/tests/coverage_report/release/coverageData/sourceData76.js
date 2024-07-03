var sourceData76 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB/mhkit/river/IO/delft_3d/delft_3d_get_all_data_points.m","RawFileContents":["function result = delft_3d_get_all_data_points(delft_3d_py_object, variable, time_index)","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%","% Returns all data points for a specific variable at a given time index from a Delft 3D netCDF","% object.","%","% Parameters","% ------------","%    delft_3d_py_object: py.netCDF4._netCDF4.Dataset","%       A netCDF python object.","%    variable: string","%       Variable in the Delft3D output","%    time_index: int","%       Index of the time variable. Points will be selected at this time index","%","% Returns","% ---------","%    result: struct","%        A struct containing the data points for the specified variable at the given time index from","%        the Delft 3D netCDF object.","%        - result.time: An array of time values","%        - result.df: A Python dataframe the output data in the original format","%        - result.*: Additional data points at the selected variable and time_index","%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","","    if ~isa(delft_3d_py_object, 'py.netCDF4._netCDF4.Dataset')","        error('MATLAB:delft_3d_get_all_data_points:InvalidInput', 'Input must be a py.netCDF4._netCDF4.Dataset object.');","    end","","    if ~isstring(variable)","        disp(class(variable));","        error('MATLAB:delft_3d_get_all_data_points:InvalidInput', 'Variable must be a string.');","    end","","    if ~isnumeric(time_index) || ~isscalar(time_index) || time_index < 1 || mod(time_index, 1) ~= 0","        error('MATLAB:delft_3d_get_all_data_points:InvalidInput', 'Time index must be a positive integer scalar.');","    end","","    % Call Python function to get all data points","    python_result = py.mhkit.river.io.d3d.get_all_data_points(delft_3d_py_object, variable, pyargs('time_index', py.int(time_index)));","","    % Convert Python dataframe to struct","    df = python_result;","    result = convert_numeric_dataframe_to_struct(df);","    result.df = df;","end",""],"CoverageDisplayDataPerLine":{"Function":{"LineNumber":1,"Hits":4,"StartColumnNumbers":0,"EndColumnNumbers":88,"ContinuedLine":false},"Statement":[{"LineNumber":28,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":62,"ContinuedLine":false},{"LineNumber":29,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":121,"ContinuedLine":false},{"LineNumber":32,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":26,"ContinuedLine":false},{"LineNumber":33,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":30,"ContinuedLine":false},{"LineNumber":34,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":96,"ContinuedLine":false},{"LineNumber":37,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":99,"ContinuedLine":false},{"LineNumber":38,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":115,"ContinuedLine":false},{"LineNumber":42,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":134,"ContinuedLine":false},{"LineNumber":45,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":23,"ContinuedLine":false},{"LineNumber":46,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":53,"ContinuedLine":false},{"LineNumber":47,"Hits":4,"StartColumnNumbers":4,"EndColumnNumbers":19,"ContinuedLine":false}]}}