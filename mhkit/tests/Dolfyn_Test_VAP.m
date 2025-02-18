classdef Dolfyn_Test_VAP < matlab.unittest.TestCase
    methods(Test)
        function test_basic_velocity_magnitude(testCase)
            % Test basic velocity magnitude calculation with simple values
            ds.vel.data = zeros(2,2,1,4);  % time, range, dir, components
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.data(:,:,:,1) = 3 * ones(2,2,1); % East component
            ds.vel.data(:,:,:,2) = 4 * ones(2,2,1); % North component
            ds.vel.dims = {'time', 'range', 'dir'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "earth";

            result = calculate_horizontal_speed_and_direction(ds);

            % Magnitude should be 5 (3-4-5 triangle)
            expected_magnitude = 5 * ones(2,2,1);
            testCase.verifyEqual(result.U_mag.data, expected_magnitude, 'AbsTol', 1e-10);
            testCase.verifyEqual(result.U_mag.units, "m s-1");
        end

        function test_different_dir_index(testCase)
            % Test with dir in different dimension positions
            ds.vel.data = zeros(2,4,2);  % time, dir, range
            ds.vel.dims = {'time', 'dir', 'range'};
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "earth";

            % Set velocity components
            idx_east = repmat({':'}, 1, ndims(ds.vel.data));
            idx_north = repmat({':'}, 1, ndims(ds.vel.data));
            idx_east{2} = 1;  % E component in dir dimension
            idx_north{2} = 2;  % N component in dir dimension

            ds.vel.data(idx_east{:}) = 3;  % East component
            ds.vel.data(idx_north{:}) = 4;  % North component

            result = calculate_horizontal_speed_and_direction(ds);

            % Verify dimensions are correct (dir dimension removed)
            testCase.verifyEqual(result.U_mag.dims, {'time', 'range'});
            % Verify magnitude calculation
            expected_magnitude = 5 * ones(2,2);  % 3-4-5 triangle
            testCase.verifyEqual(result.U_mag.data, expected_magnitude, 'AbsTol', 1e-10);
        end

        function test_squeezed_dimensions(testCase)
            % Test with singleton dimensions that need squeezing
            ds.vel.data = zeros(2,1,4,1,2);  % time, singleton, dir, singleton, range
            ds.vel.dims = {'time', 'dir', 'range'};
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "earth";

            % Set velocity components
            idx_east = repmat({':'}, 1, ndims(ds.vel.data));
            idx_north = repmat({':'}, 1, ndims(ds.vel.data));
            idx_east{3} = 1;  % E component in dir dimension
            idx_north{3} = 2;  % N component in dir dimension

            ds.vel.data(idx_east{:}) = 3;  % East component
            ds.vel.data(idx_north{:}) = 4;  % North component

            result = calculate_horizontal_speed_and_direction(ds);

            % Verify dimensions are correct after squeezing
            testCase.verifyEqual(result.U_mag.dims, {'time', 'range'});
            % Verify magnitude calculation
            expected_magnitude = 5 * ones(2,2);  % 3-4-5 triangle
            testCase.verifyEqual(result.U_mag.data, expected_magnitude, 'AbsTol', 1e-10);
        end

        function test_zero_velocity(testCase)
            % Test behavior with zero velocity
            ds.vel.data = zeros(2,2,4);
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.coord_sys = "earth";
            ds.vel.dims = {'time', 'range', 'dir'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();

            result = calculate_horizontal_speed_and_direction(ds);

            % Verify zero magnitude
            testCase.verifyEqual(sum(result.U_mag.data, 'all'), 0);
            % Direction should be NaN for zero velocity
            testCase.verifyTrue(isnan(sum(result.U_dir.data, 'all')));
        end

        function test_cardinal_directions(testCase)
            % Create base dataset
            ds.vel.data = zeros(2,2,1,4);
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.dims = {'time', 'range', 'dir'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "earth";

            % Test East direction (0 degrees mathematical, 90 degrees compass)
            ds.vel.data(:,:,:,1) = 1;  % East = 1
            ds.vel.data(:,:,:,2) = 0;  % North = 0
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 90 * ones(2,2,1), 'AbsTol', 1e-10);

            % Test North direction (90 degrees mathematical, 0 degrees compass)
            ds.vel.data(:,:,:,1) = 0;  % East = 0
            ds.vel.data(:,:,:,2) = 1;  % North = 1
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 0 * ones(2,2,1), 'AbsTol', 1e-10);

            % Test West direction (180 degrees mathematical, 270 degrees compass)
            ds.vel.data(:,:,:,1) = -1;  % East = -1
            ds.vel.data(:,:,:,2) = 0;   % North = 0
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 270 * ones(2,2,1), 'AbsTol', 1e-10);

            % Test South direction (270 degrees mathematical, 180 degrees compass)
            ds.vel.data(:,:,:,1) = 0;   % East = 0
            ds.vel.data(:,:,:,2) = -1;  % North = -1
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 180 * ones(2,2,1), 'AbsTol', 1e-10);
        end

        function test_intercardinal_directions(testCase)
            % Create base dataset
            ds.vel.data = zeros(2,2,1,4);
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.dims = {'time', 'range', 'dir'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "earth";

            % Test Northeast (45 degrees mathematical, 45 degrees compass)
            ds.vel.data(:,:,:,1) = 1;  % East = 1
            ds.vel.data(:,:,:,2) = 1;  % North = 1
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 45 * ones(2,2,1), 'AbsTol', 1e-10);

            % Test Northwest (135 degrees mathematical, 315 degrees compass)
            ds.vel.data(:,:,:,1) = -1;  % East = -1
            ds.vel.data(:,:,:,2) = 1;   % North = 1
            result = calculate_horizontal_speed_and_direction(ds);
            testCase.verifyEqual(result.U_dir.data, 315 * ones(2,2,1), 'AbsTol', 1e-10);
        end

        function test_non_earth_coordinate_system(testCase)
            % Create base dataset
            ds.vel.data = zeros(2,2,1,4);
            ds.vel.coords.dir = {'E', 'N', 'U1', 'U2'};
            ds.vel.dims = {'time', 'range', 'dir'};
            ds.vel.coords.time = 1:2;
            ds.vel.coords.range = 1:2;
            ds.attrs = struct();
            ds.coord_sys = "principal";  % Change to non-earth coordinate system

            result = calculate_horizontal_speed_and_direction(ds);

            % Verify that U_dir field is present
            testCase.verifyTrue(isfield(result, 'U_dir'));
            testCase.verifyEqual(result.U_dir.units, "degrees clockwise from E");

            % Verify that magnitude is still calculated
            testCase.verifyTrue(isfield(result, 'U_mag'));
        end

        function test_error_conditions(testCase)
            % Test missing vel field
            ds = struct('attrs', struct(), 'coord_sys', "earth");
            testCase.verifyError(@() calculate_horizontal_speed_and_direction(ds), ...
                'MATLAB:assert:failed');

            % Test missing coords.dir
            ds.vel.data = zeros(2,2,1,4);
            ds.vel.coords = struct();
            ds.vel.dims = {'time', 'range', 'dir'};
            testCase.verifyError(@() calculate_horizontal_speed_and_direction(ds), ...
                'MATLAB:assert:failed');

            % Test missing E component
            ds.vel.data = zeros(2,2,1,4);
            ds.vel.coords.dir = {'N', 'U1', 'U2'};  % Remove E
            ds.vel.dims = {'time', 'range', 'dir'};
            testCase.verifyError(@() calculate_horizontal_speed_and_direction(ds), ...
                'MATLAB:assert:failed');
        end

    end
end
