classdef River_TestIO_Delft < matlab.unittest.TestCase

    methods(Test)

        function testReadD3DFile(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);
            testCase.assertEqual(class(data), 'py.netCDF4._netCDF4.Dataset');
        end

        function testGetKeys(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);
            data_struct = delft_3d_get_keys(data);
            num_fieldnames = length(fieldnames(data_struct));
            testCase.assertEqual(data_struct.s1, 'water level');
            testCase.assertEqual(num_fieldnames, 52);
        end

    end

end
