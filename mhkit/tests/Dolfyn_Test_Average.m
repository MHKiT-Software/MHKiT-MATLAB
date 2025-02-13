classdef Dolfyn_Test_Average < matlab.unittest.TestCase

    % Define properties to store test data
    properties
        time_vec
        range_vec
        dir_vec
        ds
    end

    methods(TestMethodSetup)
        function create_sample_adcp_data(test_case)
            % Create a typical ADCP dataset structure
            test_case.time_vec = (1:24)';  % 24 hours
            test_case.range_vec = (1:28)';  % 28 depth bins
            test_case.dir_vec = (0:90:270)';  % 4 directions (typically ADCP beams)

            % Basic dataset structure
            test_case.ds.time = test_case.time_vec;
            test_case.ds.coords.time = test_case.time_vec;
            test_case.ds.range = test_case.range_vec;
            test_case.ds.coords.range = test_case.range_vec;
            test_case.ds.dir = test_case.dir_vec;
            test_case.ds.coords.dir = test_case.dir_vec;

            % Create velocity data with three dimensions [time x range x direction]
            [t, r, d] = ndgrid(1:24, 1:28, 1:4);
            test_case.ds.velocity.dims = {'time', 'range', 'dir'};
            test_case.ds.velocity.data = sin(t/12) .* cos(r/5) .* sin(d/4);
            test_case.ds.velocity.coords.time = test_case.time_vec;
            test_case.ds.velocity.coords.range = test_case.range_vec;
            test_case.ds.velocity.coords.dir = test_case.dir_vec;

            % Add echo intensity with same dimensions
            test_case.ds.echo_intensity.dims = {'time', 'range', 'dir'};
            test_case.ds.echo_intensity.data = 50 + 10*exp(-r/5) .* (1 + 0.2*sin(t/6));
            test_case.ds.echo_intensity.coords.time = test_case.time_vec;
            test_case.ds.echo_intensity.coords.range = test_case.range_vec;
            test_case.ds.echo_intensity.coords.dir = test_case.dir_vec;

            % Add temperature with only time dimension
            test_case.ds.temperature.dims = {'time'};
            test_case.ds.temperature.data = 20 + 2*sin(test_case.time_vec/12);
            test_case.ds.temperature.coords.time = test_case.time_vec;
        end
    end

    methods(Test)
        function test_time_averaging(test_case)
            % Test averaging over time dimension
            result = average_by_dimension(test_case.ds, 4, 'time');

            % Verify dimensions
            test_case.verifySize(result.velocity.data, [6 28 4]);
            test_case.verifySize(result.echo_intensity.data, [6 28 4]);
            test_case.verifySize(result.temperature.data, [6 1]);

            % Verify coordinate preservation
            test_case.verifyEqual(result.coords.time, test_case.time_vec(1:4:24));
            test_case.verifyEqual(result.coords.range, test_case.range_vec);
            test_case.verifyEqual(result.coords.dir, test_case.dir_vec);
        end

        function test_range_averaging(test_case)
            % Test averaging over range dimension
            result = average_by_dimension(test_case.ds, 2, 'range');

            % Verify dimensions
            test_case.verifySize(result.velocity.data, [24 14 4]);
            test_case.verifySize(result.echo_intensity.data, [24 14 4]);

            % Verify coordinate preservation
            test_case.verifyEqual(result.coords.time, test_case.time_vec);
            test_case.verifyEqual(result.coords.range, test_case.range_vec(1:2:28));
            test_case.verifyEqual(result.coords.dir, test_case.dir_vec);
        end

        function test_direction_averaging(test_case)
            % Test averaging over direction dimension
            result = average_by_dimension(test_case.ds, 2, 'dir');

            % Verify dimensions
            test_case.verifySize(result.velocity.data, [24 28 2]);
            test_case.verifySize(result.echo_intensity.data, [24 28 2]);

            % Verify coordinate preservation
            test_case.verifyEqual(result.coords.time, test_case.time_vec);
            test_case.verifyEqual(result.coords.range, test_case.range_vec);
            test_case.verifyEqual(result.coords.dir, test_case.dir_vec(1:2:4));
        end

        function test_nan_handling(test_case)
            % Modify velocity data to include NaNs
            ds_with_nans = test_case.ds;
            ds_with_nans.velocity.data(1:4:end, :, :) = NaN;

            result = average_by_dimension(ds_with_nans, 4, 'time');

            % Verify NaN handling
            test_case.verifyFalse(any(isnan(result.velocity.data(:))));
        end

        function test_single_variable_averaging(test_case)
            % Test averaging of temperature (single dimension)
            result = average_by_dimension(test_case.ds, 6, 'time');

            % Verify temperature averaging
            expected_size = [4 1];  % 24/6 = 4 time points
            test_case.verifySize(result.temperature.data, expected_size);

            % Verify coordinate preservation
            test_case.verifyEqual(result.temperature.coords.time, test_case.time_vec(1:6:24));
        end

        function test_uneven_samples(test_case)
            % Test with number of samples that doesn't divide evenly
            result = average_by_dimension(test_case.ds, 5, 'time');
            
            % Should use floor(24/5) = 4 complete bins
            % Other dimensions should remain unchanged (28 range bins and 4 directions)
            expected_time_points = 4;
            test_case.verifySize(result.velocity.data, [expected_time_points 28 4]);
        end

        function test_dimension_validation(test_case)
            % Test with invalid dimension name
            test_case.verifyError(@() average_by_dimension(test_case.ds, 2, 'invalid_dim'), ...
                'MATLAB:error');
        end

        function test_boundary_conditions(test_case)
            % Test with minimal dataset
            minimal_ds = struct();
            minimal_ds.time = [1; 2];
            minimal_ds.coords.time = [1; 2];
            minimal_ds.velocity.dims = {'time'};
            minimal_ds.velocity.data = [1; 2];
            minimal_ds.velocity.coords.time = [1; 2];

            result = average_by_dimension(minimal_ds, 2);
            test_case.verifySize(result.velocity.data, [1 1]);
        end
    end
end
