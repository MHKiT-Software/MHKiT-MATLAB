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

        function testGetAllTime(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            time = delft_3d_get_all_time(data);
            testCase.assertEqual(time(1, 5), 240);

            num_time = length(time);
            testCase.assertEqual(num_time, 5);
        end

        function testConvertTime(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            seconds_run = 62;
            expected_index = 2;
            closest_time_index = delft_3d_convert_time(data, seconds_run);
            testCase.assertEqual(closest_time_index, expected_index);

            seconds_run = 111;
            expected_index = 3;
            closest_time_index = delft_3d_convert_time(data, seconds_run);
            testCase.assertEqual(closest_time_index, expected_index);

            seconds_run = 1000;
            expected_index = 5;
            closest_time_index = delft_3d_convert_time(data, seconds_run);
            testCase.assertEqual(closest_time_index, expected_index);
        end

        function testGetAllDataPoints(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);
            run_timestamps = delft_3d_get_all_time(data);

            x_velocity_key = "ucx";
            seconds_into_simulation = run_timestamps(1, 4);
            seconds_into_simulation_index = delft_3d_convert_time(data, seconds_into_simulation);

            x_velocity_point_data = delft_3d_get_all_data_points(data, x_velocity_key, seconds_into_simulation_index);

            % Fieldnames
            expected_fieldnames = 7;
            num_fieldnames = length(fieldnames(x_velocity_point_data));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 2;
            expected_time = 240;
            expected_y = 1.125;

            testCase.assertEqual(x_velocity_point_data.time(1, expected_index), expected_time);
            testCase.assertEqual(x_velocity_point_data.y(1, expected_index), expected_y);

        end

    end

end
