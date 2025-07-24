classdef Dolfyn_TestIO < matlab.unittest.TestCase

    properties
        nens = 100;
    end

    methods (Test)

        % ADP Test Cases
        function test_io_rdi(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_test01.000');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_7f79(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_7f79.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_7f79.000');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_withBT(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/RDI_withBT.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/RDI_withBT.000','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_vm(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vmdas01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vmdas01.ENX', 'nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_wr1(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/winriver01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/winriver01.PD0');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_rdi_wr2(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/winriver02.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/winriver02.PD0');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        % Norteck Test Cases
        function test_io_norteck(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr','userdata',false,'nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_norteck_ud(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_ud.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_norteck_h(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/H-AWAC_test01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/H-AWAC_test01.wpr');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        % Signature Test Cases
        function test_io_signature(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/BenchFile01.nc');
            try
                delete '../../examples/data/dolfyn/BenchFile01.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/BenchFile01.ad2cp','nens',testCase.nens');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_i(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU.nc');
            try
                delete '../../examples/data/dolfyn/Sig1000_IMU.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig1000_IMU.ad2cp','userdata',false,'nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_i_ud(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_ud.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig1000_IMU.ad2cp','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_ieb(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/VelEchoBT01.nc');
            try
                delete '../../examples/data/dolfyn/VelEchoBT01.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/VelEchoBT01.ad2cp','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_ie(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo.nc');
            try
                delete '../../examples/data/dolfyn/Sig500_Echo.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig500_Echo.ad2cp','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_tide(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig1000_tidal.nc');
            try
                delete '../../examples/data/dolfyn/Sig1000_tidal.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig1000_tidal.ad2cp','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_skip(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig_SkippedPings01.nc');
            try
                delete '../../examples/data/dolfyn/Sig_SkippedPings01.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig_SkippedPings01.ad2cp');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_badt(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig1000_BadTime01.nc');
            try
                delete '../../examples/data/dolfyn/Sig1000_BadTime01.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/control/Sig1000_BadTime01.ad2cp');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_sig_leiw(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/Sig500_last_ensemble_is_whole.nc');
            try
                delete '../../examples/data/dolfyn/Sig500_last_ensemble_is_whole.ad2cp.index';
            catch
            end
            ds_read  = dolfyn_read('../../examples/data/dolfyn/Sig500_last_ensemble_is_whole.ad2cp');
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        % ADV Test Cases
        function test_io_adv(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data01.VEC','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_imu(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC','userdata',false,'nens',testCase.nens);
            warning('on','all')
            % These values are not correct for this data but I'm adding
            % them for test purposes only.
            ds_read = set_inst2head_rotmat(ds_read, eye(3));
            ds_read.attrs.inst2head_vec = [-1.0; 0.5; 0.2];
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_imu_userdata(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01-json.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        function test_io_adv_burst(testCase)
            warning('off','all')
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/control/burst_mode01.nc');
            ds_read  = dolfyn_read('../../examples/data/dolfyn/burst_mode01.VEC','nens',testCase.nens);
            warning('on','all')
            Obj.diff = Dolfyn_TestIO.compare_structures(...
                ds_read, ds_cntrl);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end

        % NetCDF Readers
        function test_io_netcdf_read_basic(testCase)
            this_ds = read_netcdf('../../examples/data/dolfyn/control/burst_mode01.nc');

            n_fieldnames = length(fieldnames(this_ds));
            expected_n_fieldnames = 19;

            testCase.assertEqual(n_fieldnames, expected_n_fieldnames);
        end

        function test_io_netcdf_read_data_attrs(testCase)
            this_ds = read_netcdf('../../examples/data/dolfyn/control/burst_mode01.nc');

            expected_pressure_unit = 'dbar';
            % TODO: Find a file with long names
            expected_pressure_long_name = 'Pressure';
            expected_temperature_unit = 'deg C';
            % TODO: Find a file with long names
            expected_temperature_long_name = 'Degrees Celsius';

            testCase.assertEqual(this_ds.pressure.units, expected_pressure_unit);
            testCase.assertEqual(this_ds.temp.units, expected_temperature_unit);
        end
    end

    methods (Static)

        function diff = compare_structures(ds_read, ds_cntrl)
            %%%%%%%%%%%%%%%%%%%%
            %     Compare the data between two structures and determine
            %     if it is within the tolerance atol.
            %
            % Parameters
            % ------------
            %     ds_read: structure
            %         Structure from the binary instrument data
            %
            %     ds_cntrl: structure
            %         Control structure read from python generated NetCDF
            %
            % Returns
            % ---------
            %     diff: float
            %         difference between the data in the two structures
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oldFmt = get(0,'Format');
            format long
            diff = 0.0;
            exclude = {'coords', 'attrs', 'time', 'complex_vars', ...
                'filehead_config'};
            % Check coords first
            fields = fieldnames(ds_read.coords);
            for qq = 1:numel(fields)
                field = fields{qq};
                if ~any(contains(field, exclude))
                    if iscell(ds_cntrl.coords.(field))
                        for kk = 1:length(ds_cntrl.coords.(field))
                            diff = diff + double(~strcmpi(...
                                ds_cntrl.coords.(field)(kk),...
                                ds_read.coords.(field)(kk)));
                        end
                    else
                        % its numeric because time got excluded
                        if size(ds_read.coords.(field)) ~= ...
                                size(ds_cntrl.coords.(field))
                            diff = (diff +  abs(...
                                sum(abs(double(ds_cntrl.coords.(field)) -...
                                double(ds_read.coords.(field)'))))/2.0);
                        else
                            diff = diff + abs(...
                                sum(abs(double(ds_cntrl.coords.(field)) -...
                                double(ds_read.coords.(field))))/2.0);
                        end
                    end
                end
            end

            % Check Attributes
            fields = fieldnames(ds_cntrl.attrs);
            for qq = 1:numel(fields)
                field = fields{qq};
                if ~any(contains(field, exclude))
                    if iscell(ds_cntrl.attrs.(field))
                        for kk = 1:numel(ds_cntrl.attrs.(field))
                            chk_nm = ds_cntrl.attrs.(field){kk};
                            diff = diff + ...
                                double(~any(strcmpi(chk_nm,...
                                ds_read.attrs.(field))));
                        end
                    elseif ischar(ds_cntrl.attrs.(field))
                        if strcmpi(ds_cntrl.attrs.(field),'yes')
                            diff = diff + double(~ds_read.attrs.(field));
                        elseif strcmpi(ds_cntrl.attrs.(field),'no')
                            diff = diff + double(ds_read.attrs.(field));
                        else
                            diff = diff + ...
                                double(~strcmpi(ds_cntrl.attrs.(field),...
                                ds_read.attrs.(field)));
                        end
                    elseif islogical(ds_read.attrs.(field))
                        double(ds_cntrl.attrs.(field) ~= ...
                                ds_read.attrs.(field));
                    elseif isnumeric(ds_cntrl.attrs.(field))
                        if contains(class(ds_cntrl.attrs.(field)),'int')
                            diff = diff + ...
                                sum(double(ds_cntrl.attrs.(field) ~=...
                                ds_read.attrs.(field)'));
                        else
                            if all(size(ds_cntrl.attrs.(field)) ==...
                                    size(ds_read.attrs.(field)))
                                atr_vl = ds_read.attrs.(field);
                            else
                                atr_vl = ds_read.attrs.(field)';
                            end
                            diff = diff + sum(abs(ds_cntrl.attrs.(field)...
                                - atr_vl), ...
                                1:numel(size(ds_cntrl.attrs.(field))))...
                                /length(ds_cntrl.attrs.(field));
                        end
                    end
                end
            end

            % Now check the remaining fields
            fields = fieldnames(ds_cntrl);
            for qq = 1:numel(fields)
                field = fields{qq};
                if ~any(contains(field, exclude))
                    cls = class(ds_cntrl.(field));
                    if strcmp(cls,'struct')
                        % Loop through structure
                        % Data
                        tmp1 = isnan(ds_cntrl.(field).data);
                        dt1 = double(ds_cntrl.(field).data);
                        dt1(tmp1) = 0.0;
                        tmp2 = isnan(ds_read.(field).data);
                        dt2 = double(ds_read.(field).data);
                        dt2(tmp2) = 0.0;

                        diff = diff + abs(sum(abs(dt1 - dt2),...
                                1:numel(size(dt1)))/length(dt1));
                        % Dims
                        for kk = 1:length(ds_cntrl.(field).dims)
                            diff = diff + double(~strcmpi(...
                                ds_cntrl.(field).dims{kk},...
                                ds_read.(field).dims{kk}));
                        end
                        % Coords (we already checked coords so just check
                        % names)
                        cntrl_coord_names =fieldnames(ds_cntrl.(field).coords);
                        read_coord_names = fieldnames(ds_read.(field).coords);
                        for kk = 1:numel(cntrl_coord_names)
                            chk_nm = cntrl_coord_names{kk};
                            diff = diff + ...
                                double(~any(strcmpi(chk_nm,read_coord_names)));
                        end
                        % Units if they exist
                        if isfield(ds_cntrl.(field),'units') && ...
                            ~strcmpi(ds_cntrl.(field).units, "none")
                            diff = diff + ...
                                double(~strcmpi(ds_cntrl.(field).units,...
                                ds_read.(field).units));
                        end
                    else
                        diff = diff + ...
                                double(~strcmpi(ds_cntrl.(field),...
                                ds_read.(field)));
                    end
                end
            end
            %fprintf('Final Diff = %f\n',diff)
            format(oldFmt);
        end

    end

end
