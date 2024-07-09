function df = convert_numeric_struct_to_dataframe(input_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Converts a MATLAB struct to a numeric pandas DataFrame.
%
% Parameters
% ------------
%    input_struct: struct
%       A MATLAB struct containing numeric arrays for each field.
%
% Returns
% ---------
%    df: DataFrame
%        A numeric pandas DataFrame containing the data from the struct.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Guard to check type of input_struct
   if ~isstruct(input_struct)
       error('MATLAB:convert_numeric_struct_to_dataframe:InvalidInput', 'input_struct must be a MATLAB struct.');
   end

   % Getting field names from the struct
   fields = fieldnames(input_struct);

   % Guard to check that each field is a numeric array
   for i = 1:numel(fields)
       field_data = input_struct.(fields{i});
       if ~isnumeric(field_data)
           error('MATLAB:convert_numeric_struct_to_dataframe:NonNumericValues', 'All values in the struct must be numeric arrays.');
       end
   end

   % Guard to check that each field value has the same length
   field_lengths = cellfun(@length, struct2cell(input_struct));
   if ~all(field_lengths == field_lengths(1))
       error('MATLAB:convert_numeric_struct_to_dataframe:FieldLengthMismatch', 'All field values must have the same length.');
   end

   % Initialize a struct to store Python lists for each field
   py_list_struct = struct();

   % Convert each array to a Python list
   for i = 1:numel(fields)
       % Access data for this field
       this_data = input_struct.(fields{i});
       % Convert the numerical array to a Python list directly
       py_list_struct.(fields{i}) = py.list(double(this_data));
   end

   % Create the DataFrame using Python
   % Convert the struct to a Python dictionary and then to a DataFrame
   df = py.pandas.DataFrame(py.dict(py_list_struct));
end
