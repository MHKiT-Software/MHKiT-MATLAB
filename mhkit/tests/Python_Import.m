% Unit tests for importing Python modules in MATLAB
classdef Python_Import < matlab.unittest.TestCase

    methods (Test)

        % Test importing 'mhkit' Python module
        function test_mhkit_import(testCase)
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

        % Test importing 'pandas' Python module
        function test_pandas_import(testCase)
            py.importlib.import_module('pandas');
            assertTrue(testCase, true);
        end

        % Test importing 'h5py' Python module
        function test_h5py_import(testCase)
            py.importlib.import_module('h5py');
            assertTrue(testCase, true);
        end

        % Verify output of the 'circular' function from 'mhkit'
        function verify_mhkit_output(testCase)
            py.importlib.import_module('mhkit');
            [x, y] = circular(30);

            expected_x = 30;
            expected_y = 706.8583;

            y_variance = abs(y - expected_y);

            assertEqual(testCase, x, expected_x);
            assertLessThan(testCase, y_variance, 0.01);
        end

        % Verify output of the 'circular' function from 'mhkit'
        function verify_mhkit_output_without_importing_module(testCase)
            % py.importlib.import_module('mhkit');
            [x, y] = circular(30);

            expected_x = 30;
            expected_y = 706.8583;

            y_variance = abs(y - expected_y);

            assertEqual(testCase, x, expected_x);
            assertLessThan(testCase, y_variance, 0.01);
        end

    end

end

