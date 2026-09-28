classdef Utils_TestMhkitUtils < matlab.unittest.TestCase
    % Tests for the native MATLAB mhkit_* input/output helpers

    methods (Test)

        function test_mhkit_standardize_user_input_to_column_vectors(testCase)
            % Row vector is transposed and flagged
            row = [1 2 3 4];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(row, 'function_name', 'test_fn');
            assertEqual(testCase, out, row(:));
            assertTrue(testCase, was_row);

            % Column vector passes through unchanged, not flagged
            col = [1; 2; 3; 4];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(col, 'function_name', 'test_fn');
            assertEqual(testCase, out, col);
            assertFalse(testCase, was_row);

            % Scalar passes through unchanged, not flagged
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(5, 'function_name', 'test_fn');
            assertEqual(testCase, out, 5);
            assertFalse(testCase, was_row);

            % Matrix (multiple column-oriented vectors) passes through unchanged
            mat = [1 2; 3 4; 5 6];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(mat, 'function_name', 'test_fn');
            assertEqual(testCase, out, mat);
            assertFalse(testCase, was_row);

            % Empty input errors with the calling function's name in the identifier
            testCase.verifyError(@() mhkit_standardize_user_input_to_column_vectors([], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');

            % 3-D input errors with the calling function's name in the identifier
            testCase.verifyError(@() mhkit_standardize_user_input_to_column_vectors(ones(2,2,2), 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');
        end

        function test_mhkit_restore_column_vectors_to_user_input(testCase)
            col = [1; 2; 3; 4];

            % was_row = true restores to a row vector
            out = mhkit_restore_column_vectors_to_user_input(col, true);
            assertEqual(testCase, out, col.');

            % was_row = false leaves the column vector unchanged
            out = mhkit_restore_column_vectors_to_user_input(col, false);
            assertEqual(testCase, out, col);

            % Round trip through standardize + restore recovers the original orientation
            row = [1 2 3 4];
            [standardized, was_row] = mhkit_standardize_user_input_to_column_vectors(row, 'function_name', 'test_fn');
            restored = mhkit_restore_column_vectors_to_user_input(standardized, was_row);
            assertEqual(testCase, restored, row);
        end

        function test_mhkit_verify_is_column_vector(testCase)
            % Column vector and scalar pass silently, numeric or datetime/duration
            mhkit_verify_is_column_vector([1;2;3], 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(5, 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(datetime(2026,1,1) + hours(0:1)', 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(hours(0:1)', 'function_name', 'test_fn');

            % Row vector and matrix both raise a function-scoped error
            testCase.verifyError(@() mhkit_verify_is_column_vector([1 2 3], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidOutput');
            testCase.verifyError(@() mhkit_verify_is_column_vector([1 2; 3 4], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidOutput');

            % Unsupported type raises a function-scoped error
            testCase.verifyError(@() mhkit_verify_is_column_vector("not numeric", 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');
        end

        function test_mhkit_frequency_to_column_names(testCase)
            frequency = [0.1; 0.2; 0.3; 0.4];
            names = mhkit_frequency_to_column_names(frequency);
            assertEqual(testCase, names, ["f_0_1000Hz"; "f_0_2000Hz"; "f_0_3000Hz"; "f_0_4000Hz"]);

            % Every generated name must be a valid MATLAB identifier
            for i = 1:numel(names)
                assertTrue(testCase, isvarname(names(i)));
            end

            % Custom decimal precision
            names2 = mhkit_frequency_to_column_names([0.1; 0.2], 'decimals', 2);
            assertEqual(testCase, names2, ["f_0_10Hz"; "f_0_20Hz"]);

            % Frequencies that collide once formatted must error, not
            % silently produce duplicate table variable names
            testCase.verifyError(@() mhkit_frequency_to_column_names([0.101; 0.102], 'decimals', 2), ...
                'MHKiT:mhkit_frequency_to_column_names:DuplicateNames');
        end

        function test_mhkit_standardize_spectrum_input_matrix(testCase)
            frequency = [0.1; 0.2; 0.3];
            spectrum = [1 4; 2 5; 3 6];
            [out_spec, out_freq, out_time, style] = mhkit_standardize_spectrum_input(spectrum, 'test_fn', frequency);

            assertEqual(testCase, out_spec, spectrum);
            assertEqual(testCase, out_freq, frequency);
            assertTrue(testCase, isempty(out_time));
            assertEqual(testCase, style, "matrix");
        end

        function test_mhkit_standardize_spectrum_input_rejects_matrix_frequency(testCase)
            % frequency must be 1-D; a caller mixing up spectrum/frequency
            % args should get a clear error here, not a confusing one later.
            spectrum = [1 4; 2 5; 3 6];
            bad_frequency = [1 2; 3 4; 5 6];
            testCase.verifyError(@() mhkit_standardize_spectrum_input(spectrum, 'my_fn', bad_frequency), ...
                'MHKiT:my_fn:InvalidOutput');
        end

        function test_mhkit_standardize_spectrum_input_struct(testCase)
            S.frequency = [0.1; 0.2; 0.3];
            S.spectrum = [1 4; 2 5; 3 6];
            S.time = datetime(2026,1,1) + hours(0:1)';

            [out_spec, out_freq, out_time, style] = mhkit_standardize_spectrum_input(S, 'test_fn');

            assertEqual(testCase, out_spec, S.spectrum);
            assertEqual(testCase, out_freq, S.frequency);
            assertEqual(testCase, out_time, S.time);
            assertEqual(testCase, style, "struct");
        end

        function test_mhkit_standardize_spectrum_input_table(testCase)
            frequency = [0.1; 0.2; 0.3];
            % Time as rows: one row per spectrum, one column per frequency bin
            T = table([1;4], [2;5], [3;6], 'VariableNames', {'f1','f2','f3'});

            [out_spec, out_freq, out_time, style] = mhkit_standardize_spectrum_input(T, 'test_fn', frequency);

            assertEqual(testCase, out_spec, [1 4; 2 5; 3 6]);
            assertEqual(testCase, out_freq, frequency);
            assertTrue(testCase, isempty(out_time));
            assertEqual(testCase, style, "table");

            % Optional 'time' variable is extracted and excluded from spectrum
            T.time = datetime(2026,1,1) + hours(0:1)';
            [out_spec2, ~, out_time2, ~] = mhkit_standardize_spectrum_input(T, 'test_fn', frequency);
            assertEqual(testCase, out_spec2, [1 4; 2 5; 3 6]);
            assertEqual(testCase, out_time2, T.time);
        end

        function test_mhkit_standardize_spectrum_input_rejects_negative_frequency(testCase)
            spectrum = [1; 2; 3];
            bad_frequency = [-0.1; 0.2; 0.3];
            testCase.verifyError(@() mhkit_standardize_spectrum_input(spectrum, 'my_fn', bad_frequency), ...
                'MHKiT:my_fn:InvalidInput');
        end

        function test_mhkit_standardize_spectrum_input_warns_on_high_frequency(testCase)
            spectrum = [1; 2; 3];
            high_frequency = [0.1; 0.2; 150];
            testCase.verifyWarning(@() mhkit_standardize_spectrum_input(spectrum, 'my_fn', high_frequency), ...
                'MHKiT:my_fn:HighFrequency');
        end

        function test_mhkit_standardize_spectrum_input_timetable(testCase)
            frequency = [0.1; 0.2; 0.3];
            time = datetime(2026,1,1) + hours(0:1)';
            TT = timetable(time, [1;4], [2;5], [3;6], 'VariableNames', {'f1','f2','f3'});

            [out_spec, out_freq, out_time, style] = mhkit_standardize_spectrum_input(TT, 'test_fn', frequency);

            assertEqual(testCase, out_spec, [1 4; 2 5; 3 6]);
            assertEqual(testCase, out_freq, frequency);
            assertEqual(testCase, out_time, time);
            assertEqual(testCase, style, "timetable");
        end

        function test_mhkit_restore_spectrum_output(testCase)
            result = [5; 2.5];
            time = datetime(2026,1,1) + hours(0:1)';

            % matrix/struct: unchanged
            assertEqual(testCase, mhkit_restore_spectrum_output(result, "matrix", "Tp"), result);
            assertEqual(testCase, mhkit_restore_spectrum_output(result, "struct", "Tp"), result);

            % table, no time
            out_table = mhkit_restore_spectrum_output(result, "table", "Tp");
            assertTrue(testCase, istable(out_table));
            assertEqual(testCase, out_table.Tp, result);

            % table, with time
            out_table_t = mhkit_restore_spectrum_output(result, "table", "Tp", time);
            assertEqual(testCase, out_table_t.time, time);
            assertEqual(testCase, out_table_t.Tp, result);

            % timetable
            out_tt = mhkit_restore_spectrum_output(result, "timetable", "Tp", time);
            assertTrue(testCase, istimetable(out_tt));
            assertEqual(testCase, out_tt.Properties.RowTimes, time);
            assertEqual(testCase, out_tt.Tp, result);
        end

    end

end
