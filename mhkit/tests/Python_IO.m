% Unit tests for data conversion to and from Python
classdef Python_IO < matlab.unittest.TestCase

    methods (Test)

        % Test creation of a multi column dataframe
        function test_simple_dataframe_creation(testCase)
            % Create sample time series and index data
            time_series = { [1, 2, 3], [4, 5, 6] };
            index = {'2023-01-01', '2023-01-02'};

            % Convert the data to a Pandas DataFrame in Python
            result = py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(time_series, index, 1);

            % Check if the result is a Pandas DataFrame
            if (isa(result,'py.pandas.core.frame.DataFrame')==1)
                % Verify the size of the output
                assertEqual(testCase, double(result.values.ndim), 2);
                assertEqual(testCase, double(result.index.size), 2);
            else
                % Fail this test if the result is not a Pandas DataFrame
                verifyFail(testCase);
            end
        end

    end

end
