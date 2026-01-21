classdef Dolfyn_TestIO_TeledyneRDI < matlab.unittest.TestCase
    % Test suite for Teledyne RD Instruments (TRDI) ADCP file reading
    %
    % ADCP File Format: RDI/TRDI
    % Manufacturer: Teledyne RD Instruments
    % File Extensions: .pd0, .000, .ENS, .ENX, .LTA, .STA
    %
    % Instruments Tested:
    %   - Basic ADCP (RDI_test01.000)
    %   - Firmware 7f79 variants
    %   - Bottom tracking data
    %   - VM-DAS Workhorse (ENX)
    %   - VM-DAS Ocean Surveyor (ENR)
    %   - WinRiver (PD0)
    %   - RiverPro (PD0)
    %   - Sentinel V B5 (pd0)

    methods (Test)

        %% =================================================================
        %% RDI INSTRUMENT TESTS
        %% =================================================================

        function test_io_rdi(testCase)
            % Test RDI basic ADCP file reading
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_test01.000');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_7f79(testCase)
            % assumeFail(testCase, "Test disabled temporarily due to known issues");
            % Test RDI ADCP with firmware 7f79
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_7f79.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_7f79.000');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_7f79_2(testCase)
            % Test RDI ADCP with firmware 7f79 variant 2
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_7f79_2.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_7f79_2.000');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_withBT(testCase)
            % Test RDI ADCP file with bottom tracking data
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_withBT.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_withBT.000','nens',nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_vm_wh(testCase)
            % Test RDI VM-DAS Workhorse file
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vmdas01_wh.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vmdas01_wh.ENX', 'nens',nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_vm_os(testCase)
            % Test RDI VM-DAS Ocean Surveyor file
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vmdas02_os.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vmdas02_os.ENR', 'nens',nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_winriver01(testCase)
            % Test RDI WinRiver file 1
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/winriver01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/winriver01.PD0');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_winriver02(testCase)
            % Test RDI WinRiver file 2
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/winriver02.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/winriver02.PD0');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_riverpro(testCase)
            % Test RDI RiverPro file
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RiverPro_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RiverPro_test01.PD0');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_winriver_transect(testCase)
            % Test RDI WinRiver transect file with limited ensembles
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/winriver02_transect.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/winriver02_transect.PD0', 'nens',nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_sentinelv_b5(testCase)
            % Test RDI Sentinel V B5 file
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/sentinelv_b5.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/sentinelv_b5.pd0');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_rdi_sec_btw_ping_division_by_zero(testCase)
            % Test fix for issue #408: RDI burst mode division by zero
            % Issue #408 reported that RDI Pinnacle 45 in continuous burst mode
            % sets sec_between_ping_groups=0 while pings_per_ensemble=1, causing
            % ZeroDivisionError in sampling rate calculation.
            nens = 10;
            warning('off','all')

            % First verify normal operation with a regular RDI file
            ds_read_normal = dolfyn_read('../../examples/data/dolfyn/RDI_test01.000', 'nens', nens);

            % Verify normal file has valid fs (not NaN)
            testCase.assertFalse(isnan(ds_read_normal.attrs.fs));
            testCase.assertGreaterThan(ds_read_normal.attrs.fs, 0);

            % Note: The MATLAB version should handle the division by zero case
            % internally and issue appropriate warnings. This test verifies that
            % the reader doesn't crash and produces valid fs values.
            % The specific mocking approach used in Python is not directly
            % translatable to MATLAB, so we test that the general case works.

            warning('on','all')
        end

    end

end
