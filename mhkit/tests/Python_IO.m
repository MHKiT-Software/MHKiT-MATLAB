% Unit tests for data conversion to and from Python
classdef PythonIO < matlab.unittest.TestCase

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

        % Test passing a list of lists to python, then reading the result
        % in python
        % Expected output is the same as the input. The final argument to
        % timeseries_to_pandas determines whethe
        % function test_convert_timeseries_to_dataframe_without_transpose(testCase)
        %     time_series = { [1, 2, 3], [4, 5, 6] };
        %     disp(class(time_series));
        %     index = {'2023-01-01', '2023-01-02'};
        %     disp(class(index));
        %     datetime_index = datetime(index, 'InputFormat', 'yyyy-MM-dd');
        %     disp(datetime_index);
        %     disp(class(datetime_index));
        %     disp(length(datetime_index));
        %     % py.mhkit_python_utils.pandas_dataframe.convert_timeseries_to_dataframe(time_series, index)
        %     df = py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(time_series, datetime_index, 1);
        %     % disp(result.values);
        %     % disp(result.index.values);
        %     % disp(class(result.index.values));
        %
        %     if (isa(df,'py.pandas.core.frame.DataFrame')==1)
        %         % disp('Is a dataframe!')
        %         % disp(result.values);
        %         % disp(class(result.values));
        %         % details(result.values);
        %         assertEqual(testCase, double(df.values.ndim), 2);
        %         assertEqual(testCase, double(df.index.size), 2);
        %
        %         df_column_count = size(df.values, 2);
        %         data=py.list();
        %         if df_column_count > 1
        %             for i = 1:x(2)
        %                 app=py.list(result.values(:,i));
        %                 data=py.mhkit_python_utils.pandas_dataframe.lis(data,app);
        %             end
        %         elseif df_column_count == 1
        %             data=df.values;
        %         end
        %
        %         disp(df.index.values);
        %         disp(class(df.index.values));
        %         % cleaned_index = cellstr(strsplit(df.index.values));
        %         % assertEqual(testCase, index, cleaned_index);
        %         % df_index = cellstr(df.index.values);
        %
        %         % Unpack the values
        %         % disp(data);
        %         cleaned_result = double(data);
        %        %  disp(cleaned_result);
        %        %
        %        % disp(class(cleaned_result));
        %        % disp(cleaned_result(1, :));
        %        % disp(cleaned_result(2, :));
        %
        %        assertEqual(testCase, time_series{1}, cleaned_result(1, :));
        %        assertEqual(testCase, time_series{2}, cleaned_result(2, :));
        %
        %     else
        %         verifyFail(testCase);
        %     end
        %
        %     assertTrue(testCase, true);
        % end

    end

end
