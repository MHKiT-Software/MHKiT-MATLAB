classdef Python_MHKiT_Import < matlab.unittest.TestCase

    methods (Test)

        function test_mhkit_import(testCase)
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

        function verify_mhkit_output(testCase)
            [x,y] = circular(30)
            
            expected_x = 30
            expected_y = 11309.733552

            y_variance = abs(y - expected_y)

            assertEqual(testCase, x, expected_x)
            assertLessThan(testCase, y_variance, 0.0001)
        end

    end
end
