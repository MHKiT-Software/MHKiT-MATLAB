classdef Tidal_TestIO < matlab.unittest.TestCase

    methods (Test) 
        
        function test_load_noaa_data(testCase)
            file_name = '../../examples/data/tidal/s08010.json';
            data = read_noaa_json(file_name);
            
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
            assertEqual(testCase,size(data.s), [92 1]);
            assertEqual(testCase,size(data.d), [92 1]);
            assertEqual(testCase,size(data.b), [92 1]);
        end
    end
end  




        