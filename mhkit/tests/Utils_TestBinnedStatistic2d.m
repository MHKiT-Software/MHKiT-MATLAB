classdef Utils_TestBinnedStatistic2d < matlab.unittest.TestCase

    methods (Test)

        function test_binned_statistic_2d_edges_and_outside_points(testCase)
            % Half-open bins with a closed last bin, elements outside the
            % edges dropped with a warning, following MHKiT-Python convention.
            % values has y bins down the rows and x bins across the columns.
            x = [0 0.5 1 1 2 2.5 3 -1]; y = [0 0 0 1 1 2 2 0]; v = 0:7;
            bins = struct('x', struct('edges', [0 1 2 3]), 'y', struct('edges', [0 1 2]));
            fn = 'mhkit_binned_statistic_2d';
            warn_id = 'MHKiT:mhkit_binned_statistic_2d:OutsideEdges';
            c = testCase.verifyWarning(@() mhkit_binned_statistic_2d(x, y, v, 'count', bins, 'function_name', fn), warn_id);
            assertEqual(testCase, c.values, [2 1 0; 0 1 3]);
            assertEqual(testCase, c.stat, 'count');
            assertEqual(testCase, c.x_edges, [0 1 2 3]);
            assertEqual(testCase, c.y_bins, [0.5 1.5]);
            warning('off', warn_id);
            cleanup = onCleanup(@() warning('on', warn_id));
            m = mhkit_binned_statistic_2d(x, y, v, 'mean', bins, 'function_name', fn);
            assertEqual(testCase, m.values, [0.5 2 NaN; NaN 3 5]);
            f = mhkit_binned_statistic_2d(x, y, v, 'frequency', bins, 'function_name', fn);
            assertEqual(testCase, f.values, [2 1 0; 0 1 3] / 8, 'AbsTol', 1e-12);
            % In-range data is warning free
            testCase.verifyWarningFree(@() mhkit_binned_statistic_2d(x(1:7), y(1:7), v(1:7), 'count', bins, 'function_name', fn));
            % Malformed grid is rejected under the caller's name
            testCase.verifyError(@() mhkit_binned_statistic_2d(x, y, v, 'count', struct('x', 1), 'function_name', 'my_fn'), 'MHKiT:my_fn:InvalidInput');
        end

        function test_binned_statistic_2d_empty_bins(testCase)
            x = 0.5; y = 0.5; v = 2;
            bins = struct('x', struct('edges', [0 1 2]), 'y', struct('edges', [0 1]));
            fn = 'mhkit_binned_statistic_2d';
            s = mhkit_binned_statistic_2d(x, y, v, 'sum', bins, 'function_name', fn);
            assertEqual(testCase, s.values, [2 0]);
            m = mhkit_binned_statistic_2d(x, y, v, 'max', bins, 'function_name', fn);
            assertEqual(testCase, m.values, [2 NaN]);
            testCase.verifyError(@() mhkit_binned_statistic_2d(10, 10, 1, 'meen', bins, 'function_name', fn), 'MATLAB:validators:mustBeMember');
        end

        function test_binned_statistic_2d_omitnan_keeps_points_aligned(testCase)
            % A NaN at index k in any one input removes element k from all
            % three, so the surviving x, y, and values stay paired by index.
            % Element k has x = k, y = 0.5, value = 100 * k. NaN at index 3
            % in x, index 7 in y, and index 9 in values.
            n = 10;
            x = (1:n)'; y = 0.5 * ones(n, 1); v = 100 * (1:n)';
            x(3) = NaN; y(7) = NaN; v(9) = NaN;
            % one x bin per element, one y bin
            bins = struct('x', struct('edges', 1:n+1), 'y', struct('edges', [0 1]));
            fn = 'mhkit_binned_statistic_2d';
            c = mhkit_binned_statistic_2d(x, y, v, 'count', bins, 'function_name', fn, 'omitnan', true);
            s = mhkit_binned_statistic_2d(x, y, v, 'sum', bins, 'function_name', fn, 'omitnan', true);
            expected_count = ones(1, n); expected_count([3 7 9]) = 0;
            expected_sum = 100 * (1:n); expected_sum([3 7 9]) = 0;
            assertEqual(testCase, c.values, expected_count);
            assertEqual(testCase, s.values, expected_sum);
        end

        function test_binned_statistic_2d_omitnan(testCase)
            % A NaN element in values, x, or y is removed by omitnan, and
            % leaves the frequency denominator
            x = [0.5 0.5 1.5 1.5 NaN]; y = [0.5 0.5 0.5 0.5 0.5]; v = [1 NaN 3 4 5];
            bins = struct('x', struct('edges', [0 1 2]), 'y', struct('edges', [0 1]));
            fn = 'mhkit_binned_statistic_2d';
            warn_id = 'MHKiT:mhkit_binned_statistic_2d:NaNInput';
            m = testCase.verifyWarning(@() mhkit_binned_statistic_2d(x, y, v, 'mean', bins, 'function_name', fn), warn_id);
            assertEqual(testCase, m.values, [NaN 3.5]);
            m = testCase.verifyWarningFree(@() mhkit_binned_statistic_2d(x, y, v, 'mean', bins, 'function_name', fn, 'omitnan', true));
            assertEqual(testCase, m.values, [1 3.5]);
            f = mhkit_binned_statistic_2d(x, y, v, 'probability', bins, 'function_name', fn, 'omitnan', true);
            assertEqual(testCase, f.values, [1 2] / 3, 'AbsTol', 1e-12);
            % 'frequency' and 'probability' are the same statistic
            g = mhkit_binned_statistic_2d(x, y, v, 'frequency', bins, 'function_name', fn, 'omitnan', true);
            assertEqual(testCase, g.values, f.values);
        end

    end

end
