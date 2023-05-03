classdef Dolfyn_Test_Clean < matlab.unittest.TestCase       
   
    methods (Test)
        function test_GN2002(testCase)
            td = dolfyn_read('../../examples/data/dolfyn/vector_data01.VEC');
            td_imu = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC');

            %TODO make GN2002, clean_fill, fill_nan_ensemble_mean
            %mask = GN2002(td.vel, npt=20);
            %td.vel = clean_fill(td.vel, mask, method='cubic', maxgap=6);
            %td.vel_clean_1D = fill_nan_ensemble_mean(td.vel[0], mask[0],
            %fs=1, window=45);
            %td.vel_clean_2D = fill_nan_ensemble_mean(td.vel, mask, fs=1,
            %window=45);

            %mask = GN2002(td_imu.vel, npt=20)
            %td_imu.vel = clean_fill(td_imu.vel, mask, method='cubic',
            %maxgap=6);

            %if make_data:
            %    save(td, 'vector_data01_GN.nc')
            %    save(td_imu, 'vector_data_imu01_GN.nc')
            %return

            %TODO compare to clean data instead of same data
            %td_cntrl = td;
            td_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/vector_data01_GN.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td, td_cntrl);
            testCase.assertLessThan(diff, 1e-6);

            %td_imu_cntrl = td_imu;
            td_imu_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/vector_data_imu01_GN.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_imu, td_imu_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_spike_thresh(testCase)
            td = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC');

            %TODO spike_thresh, clean_fill
            %mask = spike_thresh(td.vel, thresh=10);
            %td.vel = clean_fill(td.vel, mask, method='cubic', maxgap=6);

            %if make_data:
                %save(td, 'vector_data01_sclean.nc')
            %return

            %TODO compare to clean data instead of same data
            %td_cntrl = td;
            td_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/vector_data01_sclean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td, td_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_range_limit(testCase)
            td = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC');

            %TODO range_limit
            %mask = range_limit(td.vel);
            %td.vel = clean_fill(td.vel, mask, method='cubic', maxgap=6);

            %if make_data:
                %save(td, 'vector_data01_rclean.nc')
            %return

            %TODO compare to clean data instead of same data
            %td_cntrl = td;
            td_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/vector_data01_rclean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td, td_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_clean_upADCP(testCase)
            td_awac = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr');
            td_sig = dolfyn_read('../../examples/data/dolfyn/Sig1000_tidal.ad2cp');

            td_awac = find_surface_from_P(td_awac, salinity=30);
            td_awac = nan_beyond_surface(td_awac);

            set_range_offset(td_sig, 0.6);
            td_sig = find_surface_from_P(td_sig, salinity=31);
            td_sig = nan_beyond_surface(td_sig);
            td_sig = correlation_filter(td_sig, thresh=50);

            %if make_data:
            %    save(td_awac, 'AWAC_test01_clean.nc')
            %    save(td_sig, 'Sig1000_tidal_clean.nc')
            %    return
            
            %TODO compare to cleaned data instead of same data
            td_sig_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/Sig1000_tidal_clean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_sig, td_sig_cntrl);
            testCase.assertLessThan(diff, 1e-6);

            td_awac_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/AWAC_test01_clean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_awac, td_awac_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_clean_downADCP(testCase)
            ds = dolfyn_read('../../examples/data/dolfyn/Sig500_Echo.ad2cp');

            % First remove bad data
            ds.vel = val_exceeds_thresh(ds.vel, thresh = 3);
            ds.vel = fillgaps_time(ds.vel);
            ds.vel_b5 = fillgaps_time(ds.vel_b5);
            ds.vel = fillgaps_depth(ds.vel);
            ds.vel_b5 = fillgaps_depth(ds.vel_b5);

            % Then clean below seabed
            set_range_offset(ds, 0.5);
            % find_surface(ds, "thresh", 10, "nfilt", 3); 
            % find_surface doesn't exist yet? TODO
            % ds = nan_beyond_surface(ds);
            % nan_beyond_surface has error, since find_surface doesn't get
            % called yet

            % now load python comp file and compare
            %TODO compare to cleaned data instead of same data
            %ds_cntrl = ds;
            ds_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/Sig500_Echo_clean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(ds, ds_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_orient_filter(testCase)
            td_sig = dolfyn_read('../../examples/data/dolfyn/Sig1000_IMU.ad2cp');
            %td_sig = medfilt_orient(td_sig); % TODO make this function
            %td_sig = rotate2(td_sig, 'earth', inplace=True);

            td_rdi = dolfyn_read('../../examples/data/dolfyn/RDI_test01.000');
            %td_rdi = medfilt_orient(td_rdi); % TODO make this function
            %td_rdi = rotate2(td_rdi, 'earth', inplace=True);

            %if make_data:
            %    save(td_sig, 'Sig1000_IMU_ofilt.nc')
            %    save(td_rdi, 'RDI_test01_ofilt.nc')
            %return
            
            %TODO compare to cleaned data instead of same data
            %td_sig_cntrl = td_sig;
            td_sig_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/Sig1000_IMU_ofilt.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_sig, td_sig_cntrl);
            testCase.assertLessThan(diff, 1e-6);

            %td_rdi_cntrl = td_rdi;
            td_rdi_cntrl = dolfyn_read('../../examples/data/dolfyn/test_data/RDI_test01_ofilt.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_rdi, td_rdi_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end


    end

    methods (Static)

     

    end

end

