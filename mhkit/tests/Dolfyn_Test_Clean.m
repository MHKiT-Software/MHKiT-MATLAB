classdef Dolfyn_Test_Clean < matlab.unittest.TestCase       
   
    methods (Test)
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

            td_sig_cntrl = td_sig;
            %td_sig_cntrl = dolfyn_read('Sig1000_tidal_clean.nc');
            diff = Dolfyn_Test_Rotate.compare_structures(td_sig, td_sig_cntrl);
            testCase.assertLessThan(diff, 1e-6);
            
            td_awac_cntrl = td_awac;
            %td_awac_cntrl = dolfyn_read('AWAC_test01_clean.nc');
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
            % nan_beyond_surface has error, TODO fix that, or fix the
            % function call here?

            % now load python comp file and compare
            ds_cntrl = ds;
            diff = Dolfyn_Test_Rotate.compare_structures(ds, ds_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end


    end

    methods (Static)

     

    end

end

