function result_struct = convert_numeric_dataframe_to_struct(df)

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

    columns = string(py.list(df.columns.values));
    result_struct = struct();

    for i = 1:length(columns)
        this_column = columns{i};
        result_struct.(this_column) = double(df.get(this_column).values.tolist());
    end
end
