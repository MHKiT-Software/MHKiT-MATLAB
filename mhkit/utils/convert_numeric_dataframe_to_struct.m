function result_struct = convert_numeric_dataframe_to_struct(df)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Converts a numeric pandas DataFrame to a MATLAB struct.
%
% Parameters
% ------------
%    df: DataFrame
%       A numeric pandas DataFrame.
%
% Returns
% ---------
%    result_struct: struct
%        A MATLAB struct containing the data from the DataFrame.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Guard to check type of df
   if ~isa(df, 'py.pandas.core.frame.DataFrame')
       error('MATLAB:convert_numeric_dataframe_to_struct:InvalidInput', 'df must be a pandas DataFrame.');
   end

    columns = cell(df.columns.values.tolist());
    result_struct = struct();

    for i = 1:length(columns)
        this_column = char(columns{i});
        result_struct.(this_column) = double(py.numpy.array(df.get(this_column).values));
    end
end
