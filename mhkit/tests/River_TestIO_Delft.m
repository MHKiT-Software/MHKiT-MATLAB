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

        function testIndexToSeconds(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            test_index = 2;
            expected_seconds = 60;
            seconds = delft_3d_index_to_seconds(data, test_index);
            testCase.assertEqual(seconds, expected_seconds);

            test_index = 3;
            expected_seconds = 120;
            seconds = delft_3d_index_to_seconds(data, test_index);

            testCase.assertEqual(seconds, expected_seconds);
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

        function testCreatePoints(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);
            run_timestamps = delft_3d_get_all_time(data);

            x_velocity_key = "ucx";
            seconds_into_simulation = run_timestamps(1, 4);
            seconds_into_simulation_index = delft_3d_convert_time(data, seconds_into_simulation);

            x_velocity_point_data = delft_3d_get_all_data_points(data, x_velocity_key, seconds_into_simulation_index);

            xmax = min(x_velocity_point_data.x);
            xmin = max(x_velocity_point_data.x);

            ymax = min(x_velocity_point_data.y);
            ymin = max(x_velocity_point_data.y);

            waterdepth_max = min(x_velocity_point_data.waterdepth);
            waterdepth_min = max(x_velocity_point_data.waterdepth);

            x = linspace(xmin, xmax, 50); % 50 is the numpy default
            y = mean([ymax, ymin]);
            waterdepth = mean([waterdepth_max, waterdepth_min]);

            centerline_points = delft_3d_create_points(x, y, waterdepth);

            % Fieldnames
            expected_fieldnames = 4;
            num_fieldnames = length(fieldnames(centerline_points));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 2;
            expected_x = 17.5128;
            expected_y = 3;

            testCase.assertEqual(centerline_points.x(1, expected_index), expected_x, 'AbsTol',0.01);
            testCase.assertEqual(centerline_points.y(1, expected_index),expected_y, 'AbsTol',0.01);

        end

        function testInterpolateToCenterline(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);
            run_timestamps = delft_3d_get_all_time(data);

            x_velocity_key = "ucx";
            seconds_into_simulation = run_timestamps(1, 4);
            seconds_into_simulation_index = delft_3d_convert_time(data, seconds_into_simulation);

            x_velocity_point_data = delft_3d_get_all_data_points(data, x_velocity_key, seconds_into_simulation_index);

            xmax = min(x_velocity_point_data.x);
            xmin = max(x_velocity_point_data.x);

            ymax = min(x_velocity_point_data.y);
            ymin = max(x_velocity_point_data.y);

            waterdepth_max = min(x_velocity_point_data.waterdepth);
            waterdepth_min = max(x_velocity_point_data.waterdepth);

            x = linspace(xmin, xmax, 50); % 50 is the numpy default
            y = mean([ymax, ymin]);
            waterdepth = mean([waterdepth_max, waterdepth_min]);

            centerline_points = delft_3d_create_points(x, y, waterdepth);

            points = [x_velocity_point_data.x; x_velocity_point_data.y; x_velocity_point_data.waterdepth];
            values = x_velocity_point_data.(x_velocity_key);
            xi = [centerline_points.x; centerline_points.y; centerline_points.waterdepth];

            centerline_variable = delft_3d_interpolate_to_centerline(points, values, xi);

            % Value Spot Check
            expected_index = 3;
            expected_val = 0.9799;
            testCase.assertEqual(centerline_variable(1, expected_index), expected_val, 'AbsTol',0.01);
        end

        function testGetLayerData(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            x_velocity_key = "ucx";
            layer = 2;

            layer_data = delft_3d_get_layer_data(data, x_velocity_key, layer);

            % Fieldnames
            expected_fieldnames = 7;
            num_fieldnames = length(fieldnames(layer_data));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 3;
            expected_x = 0.1250;
            expected_y = 1.3750;

            testCase.assertEqual(layer_data.x(1, expected_index), expected_x, 'AbsTol',0.01);
            testCase.assertEqual(layer_data.y(1, expected_index),expected_y, 'AbsTol',0.01);
        end

        function testCalculateTurbulentIntensity(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            run_timestamps = delft_3d_get_all_time(data);
            x_velocity_key = "ucx";

            seconds_into_simulation = run_timestamps(1, 4);
            seconds_into_simulation_index = delft_3d_convert_time(data, seconds_into_simulation);

            x_velocity_point_data = delft_3d_get_all_data_points(data, x_velocity_key, seconds_into_simulation_index);

            xmax = min(x_velocity_point_data.x);
            xmin = max(x_velocity_point_data.x);

            ymax = min(x_velocity_point_data.y);
            ymin = max(x_velocity_point_data.y);

            waterdepth_max = min(x_velocity_point_data.waterdepth);
            waterdepth_min = max(x_velocity_point_data.waterdepth);

            x2 = linspace(xmin, xmax, 100);
            y_contour = linspace(ymin, ymax, 40);
            z2 = mean([waterdepth_max, waterdepth_min]);

            contour_points = delft_3d_create_points(x2, y_contour, z2);

            turbulent_intensity = delft_3d_calculate_turbulent_intensity(data, contour_points.df, true);

            % Fieldnames
            expected_fieldnames = 9;
            num_fieldnames = length(fieldnames(turbulent_intensity));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 3;
            expected_x = 17.5164;
            expected_turbulent_intensity = 4.0759;

            testCase.assertEqual(turbulent_intensity.x(1, expected_index), expected_x, 'AbsTol',0.01);
            testCase.assertEqual(turbulent_intensity.turbulent_intensity(1, expected_index), expected_turbulent_intensity, 'AbsTol',0.01);
        end

        function testCalculateVariableInterpolation(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            variables = {'turkin1', 'ucx', 'ucy', 'ucz'};

            var_interpolation = delft_3d_calculate_variable_interpolation(data, variables, "faces", "nearest");

            % Fieldnames
            expected_fieldnames = 8;
            num_fieldnames = length(fieldnames(var_interpolation));
            testCase.assertEqual(num_fieldnames, expected_fieldnames);

            % Value Spot Check
            expected_index = 4;
            expected_x = 0.6250;
            expected_ucz = 0.0046;

            testCase.assertEqual(var_interpolation.x(1, expected_index), expected_x, 'AbsTol',0.01);
            testCase.assertEqual(var_interpolation.ucz(1, expected_index), expected_ucz, 'AbsTol',0.01);
        end

        function testCalculateUnorm(testCase)
            d3d_file = "../../examples/data/river/d3d/turbineTest_map.nc";
            data = delft_3d_open_netcdf(d3d_file);

            variables = {'turkin1', 'ucx', 'ucy', 'ucz'};

            var_interpolation = delft_3d_calculate_variable_interpolation(data, variables, "faces", "nearest");

            var_interpolation.u_mag = delft_3d_calculate_unorm(var_interpolation.ucx, var_interpolation.ucy, var_interpolation.ucz);

            % Value Spot Check
            expected_index = 4;
            expected_value = 0.6942;

            testCase.assertEqual(var_interpolation.u_mag(1, expected_index), expected_value, 'AbsTol',0.01);
        end

        function testCleanupTurbulentKineticEnergy(testCase)
            input = [0 1 2 3 0.001 -0.00002 -0.01 -0.0001 5 -5];
            threshold = -0.001;
            expected_output = [0 1 2 3 0.001 0 NaN 0 5 NaN];
            output = delft_3d_cleanup_turbulent_kinetic_energy(input, threshold);

            testCase.assertEqual(output, expected_output, 'AbsTol',0.001);
        end

    end

end
