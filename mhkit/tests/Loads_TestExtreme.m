classdef Loads_TestExtreme < matlab.unittest.TestCase

    methods (Test, TestTags = {'DebuggingActions'})

        function test_mhkit_import(testCase)
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

    end

    methods (Test)

        function test_mler_coefficients(testCase)
            % create inputs and load validation data
            fpath = '../../examples/data/loads/mler_data.csv';
            validation = readtable(fpath);
            wave_freq = linspace(0,1,500);
            js = jonswap_spectrum(wave_freq,15.1,9);
            response_desired = 1;
            RAO = validation.RAO;
            % execute function
            mler = mler_coefficients(RAO, js, response_desired);
            % assertions
            assertEqual(testCase, mler.conditioned_spectrum, validation.Res_Spec, 'RelTol',0.005)
            assertEqual(testCase, mler.phase, validation.phase, 'RelTol',0.001)
        end


    end
end