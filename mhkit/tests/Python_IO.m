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

    end

end
