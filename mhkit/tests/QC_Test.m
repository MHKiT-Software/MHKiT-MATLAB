classdef QC_Test < matlab.unittest.TestCase

    methods (Test) 
        
        function test_check_corrupt(testCase)
            %  Column C has corrupt data (-999) between 7:30 and 9:30
            simple = readtable('../../examples/data/qc/simple.xlsx');
            data.values = simple.C;
            data.time = simple.Var1;
            corrupt_vals = {-999};
            for i = 1:length(data.values)
                if data.values(i) == -999
                    expect.values(i) = NaN;
                    expect.mask(i) = 0;
                else
                    expect.values(i) = data.values(i);
                    expect.mask(i) = 1;
                end
            end
            expected.values = expect.values.';
            expected.mask = int64(expect.mask.');
            results = check_corrupt(data,corrupt_vals);
            assertEqual(testCase, results.values, expected.values);
            assertEqual(testCase, results.mask, expected.mask);
        end
              
        function test_check_delta(testCase)
            % Column A has the same value (0.5) from 12:00 until 14:30
            % Column C does not follow the expected sine function from 13:00 until 16:15. 
            % The change is abrupt and gradually corrected.
            simple = readtable('../../examples/data/qc/simple_expected.xlsx');
            simple_expected = readtable('../../examples/data/qc/simple_expected.xlsx');
            
            dataA.values = simple.A;
            dataC.values = simple.C;
            dataA.time = simple.Var1;
            dataC.time = simple.Var1;
            format long
            datenumA = datenum(dataA.time);
            datenumC = datenum(dataC.time);
            bound = [-1.0, 1.0];
            window = 2*3600; % seconds
            resultsA = check_delta(dataA,bound,window);
            resultsC = check_delta(dataC,bound,window);            
            
            expectedA.values = simple_expected.A;
            expectedC.values = simple_expected.C;
            ABSTOL = 0.00000001;
            assertEqual(testCase, resultsA.values, expectedA.values, 'AbsTol', ABSTOL);
            assertEqual(testCase, resultsC.values, expectedC.values, 'AbsTol', ABSTOL);
        end

        function test_check_increment(testCase)
            % Column A has the same value (0.5) from 12:00 until 14:30
            % Column C does not follow the expected sine function from 13:00 until 16:15.
            % The change is abrupt and gradually corrected.
            simple = readtable('../../examples/data/qc/simple.xlsx');
            dataA.values = simple.A;
            dataC.values = simple.C;
            dataA.time = simple.Var1;
            dataC.time = simple.Var1;

            format long
            datenumA = datenum(dataA.time);
            datenumC = datenum(dataC.time);

            bound = [0.0001, 0.6];
%             
            expectA.values = zeros(size(dataA.values));
            for i = 1:96
                if i >= 68 && i <= 69
                    expectA.values(i) = NaN;
                    expectA.mask(i) = 0;
                elseif i >= 49 && i <= 58
                    
                    expectA.values(i) = NaN;
                    expectA.mask(i) = 0;
                else
                    expectA.values(i) = dataA.values(i);
                    expectA.mask(i) = 1;
                end   
            end
            expectC.values = zeros(size(dataC.values));
            for ii = 1:96
                if ii >= 68 && ii <= 69
                    expectC.values(ii) = NaN;
                    expectC.mask(ii) = 0;
                elseif ii >= 30 && ii <= 39
                    expectC.values(ii) = NaN;
                    expectC.mask(ii) = 0;
                elseif ii == 52
                    expectC.values(ii) = NaN;
                    expectC.mask(ii) = 0;
                else
                    
                    expectC.values(ii) = dataC.values(ii);
                    expectC.mask(ii) = 1;
                end   
            end
            expectedA.values = expectA.values;
            expectedA.mask = int64(expectA.mask.');
            expectedC.values = expectC.values;
            expectedC.mask = int64(expectC.mask.');

            resultsA = check_increment(dataA,bound);
            resultsC = check_increment(dataC,bound);
            assertEqual(testCase, resultsA.values, expectedA.values);
            assertEqual(testCase, resultsA.mask, expectedA.mask);
            assertEqual(testCase, resultsC.values, expectedC.values);
            assertEqual(testCase, resultsC.mask, expectedC.mask);
        end

        function test_check_missing(testCase)
            % Column D is missing data from 17:45 until 18:15
            simple = readtable('../../examples/data/qc/simple.xlsx');
            data.values = simple.D;
            data.time = simple.Var1;
            A = ismissing(data.values);
            for i = 1:length(A)
                if A(i) == 1
                    expect.values(i) = NaN;
                    expect.mask(i) = 0;
                else
                    expect.values(i) = data.values(i);
                    expect.mask(i) = 1;
                end
            end
            expected.values = expect.values.';
            expected.mask = int64(expect.mask.');
            results = check_missing(data);
            assertEqual(testCase, results.values, expected.values);
            assertEqual(testCase, results.mask, expected.mask);
        end

%         function test_check_outlier(testCase)
% 
%         end  
        
%         function test_check_range(testCase)
%             % Column B is below the expected lower bound of 0 at 6:30 and above the expected upper bound of 1 at 15:30
%             % Column D is occasionally below the expected lower bound of -1 around midday (2 time steps)
%             % and above the expected upper bound of 1 in the early morning and late evening (10 time steps).
%         end
        
      

        function test_check_timestamp(testCase)
%             % Missing timestamp at 5:00
%             % Duplicate timestamp 17:00
%             % Non-monotonic timestamp 19:30
            simple = readtable('../../examples/data/qc/simple.xlsx');
            simple_expected = readtable('../../examples/data/qc/simple_expected.xlsx');
            data.values = simple.A;
            data.time = simple.Var1;
            freq = 900; % seconds
            expected.values = simple_expected.A;
            expected.time = (simple_expected.Var1);
            results = check_timestamp(data,freq);
            datenum_time = datenum(expected.time) - datenum(results.time.');
            expected_datenum_time = zeros(length(datenum_time),1);
            ABSTOL = 0.00000001;
            assertEqual(testCase, results.values, expected.values, 'AbsTol', ABSTOL);
            assertEqual(testCase, datenum_time, expected_datenum_time, 'AbsTol', ABSTOL);
        end
    end
end  
