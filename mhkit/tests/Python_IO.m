% Unit tests for data conversion to and from Python
classdef Python_IO < matlab.unittest.TestCase

    methods (Test)

         % Test creation of a Pandas DataFrame with timeseries index
        function testConvertNumericDataframeToStructure(testCase)
            % Create a Pandas DataFrame in Python
            data = { [1, 2, 3], [4, 5, 6] };
            columns = {'first', 'second', 'third'};
            df = py.pandas.DataFrame(data, pyargs('columns', columns));

            % Convert the DataFrame to a struct
            output_struct = convert_numeric_dataframe_to_struct(df);

            % Fieldnames Check
            expected_fieldnames = 3;
            num_fieldnames = length(fieldnames(output_struct));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 2;
            expected_first = 4;
            expected_second = 5;

            testCase.assertEqual(output_struct.first(1, expected_index), expected_first, 'AbsTol',0.01);
            testCase.assertEqual(output_struct.second(1, expected_index), expected_second, 'AbsTol',0.01);
        end

         % Convert a numerical struct to a pd.DataFrame
        function testConvertStructToDataFrame(testCase)
            % Create a numeric structure
            ds = struct('data1', [1, 2, 3, 4, 5], 'data2', [6, 7, 8, 9, 10], 'data3', [11, 12, 13, 14, 15]);

            % Convert the DataFrame to a struct
            output_df = convert_numeric_struct_to_dataframe(ds);

            % % Value Spot Check
            expected_index = 1;
            expected_first = 2;

            data_1 = output_df.get("data1");
            actual_first = data_1.get(expected_index);

            testCase.assertEqual(actual_first, expected_first, 'AbsTol',0.01);
        end

        function testConvertStructToDataFrameLengthError(testCase)
            % Create a numeric structure
            ds = struct('data1', [1, 2, 3, 4], 'data2', [6, 7, 8, 9, 10], 'data3', [11, 12, 13, 14, 15]);

            testCase.verifyError(@() convert_numeric_struct_to_dataframe(ds), 'MATLAB:convert_numeric_struct_to_dataframe:FieldLengthMismatch');
        end

        function testConvertStructToDataFrameTypeError(testCase)
            % Create a numeric structure
            ds = struct('data1', [1, 2, 3, 4, 5], 'data2', {"a", "b", "c", "d", "e"}, 'data3', [11, 12, 13, 14, 15]);

            testCase.verifyError(@() convert_numeric_struct_to_dataframe(ds), 'MATLAB:convert_numeric_struct_to_dataframe:NonNumericValues');
        end

    end

end
