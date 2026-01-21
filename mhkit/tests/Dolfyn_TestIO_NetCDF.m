classdef Dolfyn_TestIO_NetCDF < matlab.unittest.TestCase
    % Test suite for Dolfyn NetCDF file I/O functionality
    %
    % File Format: NetCDF
    % Standard: Network Common Data Form
    % File Extensions: .nc

    methods (Test)

        %% =================================================================
        %% NETCDF I/O TESTS
        %% =================================================================

        function test_io_netcdf_read_basic(testCase)
            % Test basic NetCDF reading functionality
            this_ds = read_netcdf('../../examples/data/dolfyn/control/vector_burst_mode01.nc');

            n_fieldnames = length(fieldnames(this_ds));
            expected_n_fieldnames = 19;

            testCase.assertEqual(n_fieldnames, expected_n_fieldnames);
        end

        function test_io_netcdf_read_data_attrs(testCase)
            % Test NetCDF data attributes reading
            this_ds = read_netcdf('../../examples/data/dolfyn/control/vector_burst_mode01.nc');

            expected_pressure_unit = 'dbar';
            expected_pressure_long_name = 'Pressure';
            expected_temperature_unit = 'degree_C';
            expected_temperature_long_name = 'Degrees Celsius';

            testCase.assertEqual(this_ds.pressure.units, expected_pressure_unit);
            testCase.assertEqual(this_ds.temp.units, expected_temperature_unit);
        end

    end

end
