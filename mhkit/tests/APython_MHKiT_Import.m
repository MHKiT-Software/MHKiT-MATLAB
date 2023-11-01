classdef APython_MHKiT_Import < matlab.unittest.TestCase

    methods (Test)

        function test_pandas_import(testCase)
            py.importlib.import_module('pandas');
            assertTrue(testCase, true);
        end

        function test_mhkit_import(testCase)
            pyenv()
            getenv('PYTHONHOME')
            getenv("PATH")
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

        function verify_mhkit_output(testCase)
            [x,y] = circular(30);
            
            expected_x = 30;
            expected_y = 11309.733552;

            y_variance = abs(y - expected_y);

            assertEqual(testCase, x, expected_x)
            assertLessThan(testCase, y_variance, 0.01)
        end

        function test_h5py_import(testCase)
            py.importlib.import_module('h5py');
            assertTrue(testCase, true);
        end

    end
end
