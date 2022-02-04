classdef Tidal_TestIO < matlab.unittest.TestCase

    methods (Test) 
        
        function test_load_noaa_data(testCase)
            
            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            
            assertTrue(testCase, isfield(data,'s'));
            assertTrue(testCase, isfield(data,'d'));
            assertTrue(testCase, isfield(data,'b')); 
            assertEqual(testCase,size(data.s), [18890 1]);
            assertEqual(testCase,size(data.d), [18890 1]);
            assertEqual(testCase,size(data.b), [18890 1]);
            
        end

        function test_request_noaa_data(testCase)
            data = request_noaa_data('s08010', 'currents','20180101','20180102');
            
            assertTrue(testCase, isfield(data,'s'));
            assertTrue(testCase, isfield(data,'d'));
            assertTrue(testCase, isfield(data,'b'));
            assertEqual(testCase,size(data.s), [184 1]);
            assertEqual(testCase,size(data.d), [184 1]);
            assertEqual(testCase,size(data.b), [184 1]);
        end
    end
end  




        