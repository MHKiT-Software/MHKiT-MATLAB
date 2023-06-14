classdef Dolfyn_Test_Clean < matlab.unittest.TestCase       
   
    methods (Test)
        function test_GN2002(testCase)
            td = dolfyn_read('../../examples/data/dolfyn/vector_data01.VEC');
            %td_imu = dolfyn_read('../../examples/data/dolfyn/vector_data_imu01.VEC');       

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
            td_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/vector_data01_GN.nc');
            diff = compare_structures(td, td_cntrl);
            testCase.assertLessThan(diff, 1e-6);

            td_imu_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/vector_data_imu01_GN.nc');
            diff = compare_structures(td_imu, td_imu_cntrl);
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
            td_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/vector_data01_sclean.nc');
            diff = compare_structures(td, td_cntrl);
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
            td_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/vector_data01_rclean.nc');
            diff = compare_structures(td, td_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_clean_upADCP(testCase)
            td_awac = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr');
            %td_sig = dolfyn_read('../../examples/data/dolfyn/Sig1000_tidal.ad2cp');

            td_awac = find_surface_from_P(td_awac, salinity=30);
            td_awac = nan_beyond_surface(td_awac);

            %set_range_offset(td_sig, 0.6);
            %td_sig = find_surface_from_P(td_sig, salinity=31);
            %td_sig = nan_beyond_surface(td_sig);
            %td_sig = correlation_filter(td_sig, thresh=50);

            %if make_data:
            %    save(td_awac, 'AWAC_test01_clean.nc')
            %    save(td_sig, 'Sig1000_tidal_clean.nc')
            %    return
            
            %TODO compare to cleaned data instead of same data
            %td_sig_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/Sig1000_tidal_clean.nc');
            %diff = compare_structures(td_sig, td_sig_cntrl);
            %testCase.assertLessThan(diff, 1e-6);

            td_awac_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/AWAC_test01_clean.nc');
            %td_awac_cntrl = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr');
            td_awac.vel
            td_awac_cntrl.vel
            diff = Dolfyn_Test_Clean.compare_structures(td_awac, td_awac_cntrl);
            testCase.assertLessThan(diff, 1e-6);
        end

        function test_clean_downADCP(testCase)
            %ds = dolfyn_read('../../examples/data/dolfyn/Sig500_Echo.ad2cp');
            ds = read_netcdf('../../examples/data/dolfyn/test_data/Sig500_Echo.nc');

            % First remove bad data
            ds.vel = val_exceeds_thresh(ds.vel, thresh = 3);
            ds.vel = fillgaps_time(ds.vel);
            ds.vel_b5 = fillgaps_time(ds.vel_b5);
            ds.vel = fillgaps_depth(ds.vel);
            ds.vel_b5 = fillgaps_depth(ds.vel_b5);

            % Then clean below seabed
            set_range_offset(ds, 0.5);
            find_surface(ds, "thresh", 10, "nfilt", 3); 
            ds = nan_beyond_surface(ds);

            % now load python comp file and compare
            ds_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/Sig500_Echo_clean.nc');
            ds.depth
            ds_cntrl.depth
            diff = Dolfyn_Test_Clean.compare_structures(ds, ds_cntrl);
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
            td_sig_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/Sig1000_IMU_ofilt.nc');
            diff = compare_structures(td_sig, td_sig_cntrl);
            testCase.assertLessThan(diff, 1e-6);

            %td_rdi_cntrl = td_rdi;
            td_rdi_cntrl = read_netcdf('../../examples/data/dolfyn/test_data/RDI_test01_ofilt.nc');
            diff = compare_structures(td_rdi, td_rdi_cntrl);
            testCase.assertLessThan(diff, 1e-6);
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
                'filehead_config', 'c_sound', 'temp', 'pressure', 'mag', ...
                'accel', 'batt', 'temp_clock', 'error', 'status', ...
                'ensemble', 'heading', 'pitch','roll'};
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
                        tmp2 = isnan(ds_read.(field).data);
                        field
                        tmp = tmp1|tmp2;

                        dt1 = double(ds_cntrl.(field).data);
                        dt1(tmp) = 0.0;                        
                        dt2 = double(ds_read.(field).data);
                        dt2(tmp) = 0.0;                       
                          
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

