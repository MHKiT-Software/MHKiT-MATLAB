classdef Dolfyn_TestIO < matlab.unittest.TestCase
    
    methods (Test)

        % ADP Test Cases  
        function test_io_rdi(testCase)
            warning('off','all')
            dat_rdi = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            dat_rdi_7f79 = read_netcdf('../../examples/data/dolfyn/control/RDI_7f79.nc');
            dat_rdi_bt = read_netcdf('../../examples/data/dolfyn/control/RDI_withBT.nc');
            dat_rdi_vm = read_netcdf('../../examples/data/dolfyn/control/vmdas01.nc');
            dat_wr1 = read_netcdf('../../examples/data/dolfyn/control/winriver01.nc');
            dat_wr2 = read_netcdf('../../examples/data/dolfyn/control/winriver02.nc');
    
            nens = 100;
            td_rdi = dolfyn_read('../../examples/data/dolfyn/RDI_test01.000');
            td_7f79 = dolfyn_read('../../examples/data/dolfyn/RDI_7f79.000');
            td_rdi_bt = dolfyn_read('../../examples/data/dolfyn/RDI_withBT.000', nens=nens);
            td_vm = dolfyn_read('../../examples/data/dolfyn/vmdas01.ENX', nens=nens);
            td_wr1 = dolfyn_read('../../examples/data/dolfyn/winriver01.PD0');
            td_wr2 = dolfyn_read('../../examples/data/dolfyn/winriver02.PD0');
            warning('on','all')
            
            assert(compare_structures(td_rdi, dat_rdi) < 1e-6);
            assert(compare_structures(td_7f79, dat_rdi_7f79) < 1e-6);
            assert(compare_structures(td_rdi_bt, dat_rdi_bt) < 1e-6);
            assert(compare_structures(td_vm, dat_rdi_vm) < 1e-6);
            assert(compare_structures(td_wr1, dat_wr1) < 1e-6);
            assert(compare_structures(td_wr2, dat_wr2) < 1e-6);
        end

        function test_io_norteck(testCase)  
            warning('off','all')
            dat_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01.nc');
            dat_awac_ud = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_ud.nc');
            dat_hwac = read_netcdf('../../examples/data/dolfyn/control/H-AWAC_test01.nc');
            
            nens = 100;
            td_awac = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr', userdata=false, nens=nens);
            td_awac_ud = dolfyn_read('../../examples/data/dolfyn/AWAC_test01.wpr', nens=nens);
            td_hwac = dolfyn_read('../../examples/data/dolfyn/H-AWAC_test01.wpr');
            warning('on','all')
            
            assert(compare_structures(td_awac, dat_awac) < 1e-6);
            assert(compare_structures(td_awac_ud, dat_awac_ud) < 1e-6);
            assert(compare_structures(td_hwac, dat_hwac) < 1e-6);
        end
        
        function test_io_signature(testCase)
            warning('off','all')
            dat_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01.nc');
            dat_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU.nc');
            dat_sig_i_ud = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_ud.nc');
            dat_sig_ieb = read_netcdf('../../examples/data/dolfyn/control/VelEchoBT01.nc');
            dat_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo.nc');
            dat_sig_tide = read_netcdf('../../examples/data/dolfyn/control/Sig1000_tidal.nc');
            dat_sig_skip = read_netcdf('../../examples/data/dolfyn/control/Sig_SkippedPings01.nc');
            dat_sig_badt = read_netcdf('../../examples/data/dolfyn/control/Sig1000_BadTime01.nc');
            dat_sig5_leiw = read_netcdf('../../examples/data/dolfyn/control/Sig500_last_ensemble_is_whole.nc');
            
            nens = 100;
            td_sig = dolfyn_read('../../examples/data/dolfyn/BenchFile01.ad2cp', nens=nens);
            td_sig_i = dolfyn_read('../../examples/data/dolfyn/Sig1000_IMU.ad2cp', userdata=false, nens=nens);
            td_sig_i_ud = dolfyn_read('../../examples/data/dolfyn/Sig1000_IMU.ad2cp', nens=nens);
            td_sig_ieb = dolfyn_read('../../examples/data/dolfyn/VelEchoBT01.ad2cp', nens=nens);
            td_sig_ie = dolfyn_read('../../examples/data/dolfyn/Sig500_Echo.ad2cp', nens=nens);
            td_sig_tide = dolfyn_read('../../examples/data/dolfyn/Sig1000_tidal.ad2cp', nens=nens);
            td_sig_skip = dolfyn_read('../../examples/data/dolfyn/Sig_SkippedPings01.ad2cp');
            td_sig_badt = dolfyn_read('../../examples/data/dolfyn/control/Sig1000_BadTime01.ad2cp'); 
            td_sig5_leiw = dolfyn_read('../../examples/data/dolfyn/Sig500_last_ensemble_is_whole.ad2cp');           
            warning('on','all')
    
            assert(compare_structures(td_sig, dat_sig) < 1e-6); 		
            assert(compare_structures(td_sig_i, dat_sig_i) < 1e-6); 		
            assert(compare_structures(td_sig_i_ud, dat_sig_i_ud) < 1e-6); 	
            assert(compare_structures(td_sig_ieb, dat_sig_ieb) < 1e-6); 	
            assert(compare_structures(td_sig_ie, dat_sig_ie) < 1e-6); 		
            assert(compare_structures(td_sig_tide, dat_sig_tide) < 1e-6); 	
            assert(compare_structures(td_sig_skip, dat_sig_skip) < 1e-6); 	
            assert(compare_structures(td_sig_badt, dat_sig_badt) < 1e-6); 	
            assert(compare_structures(td_sig5_leiw, dat_sig5_leiw) < 1e-6);
        end
        
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
            oldFmt = format;
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
                            diff = diff + abs(... 
                                sum(double(round(ds_cntrl.attrs.(field),4) -...
                                round(ds_read.attrs.(field)',4))));
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
                        
                        % The difference between significant figures in pi,
                        % cos, sin, etc between numpy and matlab causes
                        % an inherent difference in the answers for
                        % orientmat and others    
                        diff = diff + fix(abs(sum(abs(dt1 - dt2),...
                                1:numel(size(dt1)))/length(dt1))*1e4)/1e4; 
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
                        if isfield(ds_cntrl.(field),'units')
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

