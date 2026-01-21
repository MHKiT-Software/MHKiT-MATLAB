classdef Dolfyn_TestIO_NortekClassic < matlab.unittest.TestCase
    % Test suite for Nortek Classic Instruments ADCP/ADV file reading
    %
    % ADCP File Format: Nortek Classic
    % Manufacturer: Nortek
    % File Extensions: .vec, .aqd, .wpr, .hdr, .sen
    %
    % Instruments Tested:
    %   - AWAC (Acoustic Wave and Current profiler)
    %   - HR-AWAC (High Resolution AWAC)
    %   - Vector ADV (Acoustic Doppler Velocimeter)
    %   - Aquadopp
    %   - Continental

    methods (Test)

        %% =================================================================
        %% NORTEK AWAC INSTRUMENT TESTS
        %% =================================================================

        function test_io_nortek_awac(testCase)
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr', 'userdata', false, 'nens', nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_nortek_awac_userdata(testCase)
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_ud.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr', 'nens', nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_nortek_hawac(testCase)
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/H-AWAC_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/H-AWAC_test01.wpr');
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        %% =================================================================
        %% NORTEK VECTOR ADV INSTRUMENT TESTS
        %% =================================================================

        function test_io_adv_vector(testCase)
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data01.VEC', 'nens', nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_vector_imu(testCase)
            % Test Nortek Vector ADV with IMU, no userdata
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC', 'userdata', false, 'nens', nens);
            warning('on','all')
            % These values are added for test purposes to match control data
            ds_read = set_inst2head_rotmat(ds_read, eye(3));
            ds_read.attrs.inst2head_vec = [-1.0; 0.5; 0.2];
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_vector_imu_userdata(testCase)
            % Test Nortek Vector ADV with IMU and userdata
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01-json.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC', 'nens', nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_vector_burst(testCase)
            % Test Nortek Vector burst mode ADV
            nens = 100;
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_burst_mode01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_burst_mode01.VEC', 'nens', nens);
            warning('on','all')
            Obj.diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

    end

end
