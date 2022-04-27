classdef Dolfyn_Test_Rotate < matlab.unittest.TestCase       
   
    methods (Test)

        % ADV Rotation Cases  
        function test_heading(testCase) 
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc'); 
            [head, pitch, roll] = orient2euler(td);
            td.pitch.data = pitch;
            td.roll.data = roll;
            td.heading.data = head;

            cd = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_head_pitch_roll.nc');
            
            Obj.diff = Dolfyn_Test_Rotate.compare_structures(...
                td, cd);
            testCase.assertLessThan(Obj.diff, 1e-6);
        end    

        function test_inst2head_rotmat(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');

            % Swap x, y, reverse z
            tr = set_inst2head_rotmat(td,[0,1,0;1,0,0;0,0,-1]);
            testCase.assertLessThan(sum(td.vel.data(:,:,1) -...
                tr.vel.data(:,:,2)), 1e-6);
            testCase.assertLessThan(sum(td.vel.data(:,:,2) -...
                tr.vel.data(:,:,1)), 1e-6);
            testCase.assertLessThan(sum(td.vel.data(:,:,3) +...
                tr.vel.data(:,:,3)), 1e-6);

            % Validation for non-symmetric rotations
            tr = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            R = squeeze(euler2orient(20, 30, 60));
            tr = set_inst2head_rotmat(tr, R);
            vel1 = tr.vel.data;
            vel2 = reshape(squeeze(td.vel.data)*R,[100,1,3]);
            testCase.assertLessThan(sum(vel1-vel2), 1e-6);
        end

        function adv_rotate_inst2earth(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            td = rotate2(td,'earth');
            tdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc'); 
            tdm = rotate2(tdm, 'earth');
            tdo = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            omat = tdo.orientmat;
            tdo = rmfield(tdo,'orientmat');
            tdo = rotate2(tdo, 'earth');
            tdo.orientmat = omat;

            cd = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2earth.nc');
            cdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_inst2earth.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdm, cdm);
            diff3 = Dolfyn_Test_Rotate.compare_structures(tdo, cd);

            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6);
            testCase.assertLessThan(diff3, 1e-6);
        end

        function test_rotate_earth2inst(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2earth.nc');
            td = rotate2(td, 'inst');
            tdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_inst2earth.nc');
            tdm = rotate2(tdm,'inst');

            cd = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            cdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');
            % The heading/pitch/roll data gets modified during rotation, 
            % so it doesn't go back to what it was.
            cdm = rmfield(cdm, {'heading', 'pitch', 'roll'});
            tdm = rmfield(tdm, {'heading', 'pitch', 'roll'});
            
            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdm, cdm);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6);
        end

        function adv_rotate_inst2beam(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            td = rotate2(td, 'beam');
            tdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');
            tdm = rotate2(tdm, 'beam');

            cd  = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2beam.nc');
            cdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_inst2beam.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdm, cdm);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6); 
        end

        function adv_rotate_beam2inst(testCase)
            td  = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2beam.nc');
            td = rotate2(td, 'inst');
            tdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_inst2beam.nc');
            tdm = rotate2(tdm, 'inst');

            cd = read_netcdf('../../examples/data/dolfyn/control/vector_data01.nc');
            cdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdm, cdm);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6); 
        end

        function test_rotate_earth2principal(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2earth.nc');
            td.attrs.principal_heading = ...
                calc_principal_heading(td.vel.data,true);
            td = rotate2(td, 'principal');
            tdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_inst2earth.nc');
            tdm.attrs.principal_heading = ...
                calc_principal_heading(tdm.vel.data,true);
            tdm = rotate2(tdm, 'principal');

            cd = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_earth2principal.nc');
            cdm = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01_rotate_earth2principal.nc');
            
            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdm, cdm);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6); 
        end

        function test_rotate_earth2principal_set_declination(testCase)
            declin = 3.875;
            td = read_netcdf('../../examples/data/dolfyn/control/vector_data01_rotate_inst2earth.nc');
            td0 = td;

            td.attrs.principal_heading = ...
                calc_principal_heading(td.vel.data,true);
            td = rotate2(td, 'principal');
            td = set_declination(td, declin);
            td = rotate2(td,'earth');

            td0 = set_declination(td0, -1);
            td0 = set_declination(td0, declin);
            td0.attrs.principal_heading = ...
                calc_principal_heading(td0.vel.data,true);
            td0 = rotate2(td0, 'earth');
            diff = Dolfyn_Test_Rotate.compare_structures(td0, td);
            testCase.assertLessThan(diff, 1e-6);
        end

        %##################################################################

        % ADP Rotation Cases  
        function adp_rotate_beam2inst(testCase)
            td_rdi = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            td_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01.nc');
            td_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU.nc');
            td_sig_ieb = read_netcdf('../../examples/data/dolfyn/control/VelEchoBT01.nc');

            td_rdi = rotate2(td_rdi, 'inst');
            td_sig = rotate2(td_sig, 'inst');
            td_sig_i = rotate2(td_sig_i, 'inst');
            td_sig_ieb = rotate2(td_sig_ieb, 'inst');

            cd_rdi = read_netcdf('../../examples/data/dolfyn/control/RDI_test01_rotate_beam2inst.nc');
            cd_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01_rotate_beam2inst.nc');
            cd_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_rotate_beam2inst.nc');
            cd_sig_ieb = read_netcdf('../../examples/data/dolfyn/control/VelEchoBT01_rotate_beam2inst.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td_rdi, cd_rdi);
            diff2 = Dolfyn_Test_Rotate.compare_structures(td_sig, cd_sig);
            diff3 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig_i, cd_sig_i);
            diff4 = Dolfyn_Test_Rotate.compare_structures(...
                td_sig_ieb, cd_sig_ieb);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6);            
            testCase.assertLessThan(diff3, 1e-6);
            testCase.assertLessThan(diff4, 1e-6);
        end

        function adp_rotate_inst2beam(testCase)
            td = read_netcdf('../../examples/data/dolfyn/control/RDI_test01_rotate_beam2inst.nc');
            td = rotate2(td, 'beam');

            td_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_earth2inst.nc');
            td_awac = rotate2(td_awac, 'beam');

            td_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01_rotate_beam2inst.nc');
            td_sig = rotate2(td_sig, 'beam');

            td_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_rotate_beam2inst.nc');
            td_sig_i = rotate2(td_sig_i, 'beam');

            td_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo_earth2inst.nc');
            td_sig_ie = rotate2(td_sig_ie, 'beam');

            cd_td = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            cd_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_inst2beam.nc');
            cd_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01.nc');
            cd_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU.nc');
            cd_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo_inst2beam.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd_td);
            diff2 = Dolfyn_Test_Rotate.compare_structures(td_awac,cd_awac);
            diff3 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig, cd_sig);
            diff4 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig_i,cd_sig_i);
            diff5 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig_ie,cd_sig_ie);
            testCase.assertLessThan(diff1, 1e-6);
            testCase.assertLessThan(diff2, 1e-6);            
            testCase.assertLessThan(diff3, 1e-6);
            testCase.assertLessThan(diff4, 1e-6);
            testCase.assertLessThan(diff5, 1e-6);
        end

        function adp_rotate_inst2earth(testCase)
            % AWAC & Sig500 are loaded in earth
            td_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01.nc');
            td_awac = rotate2(td_awac, 'inst');
            td_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo.nc');
            td_sig_ie = rotate2(td_sig_ie, 'inst');
            td_sig_o = td_sig_ie;
    
            td = read_netcdf('../../examples/data/dolfyn/control/RDI_test01.nc');
            td = rotate2(td, 'earth');
            tdwr2 = read_netcdf('../../examples/data/dolfyn/control/winriver02.nc');
            tdwr2 = rotate2(tdwr2, 'earth');
            td_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01_rotate_beam2inst.nc');
            td_sig = rotate2(td_sig, 'earth');
            td_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_rotate_beam2inst.nc');
            warning('off','all')
            td_sig_i = rotate2(td_sig_i, 'earth');
            warning('on','all')

            td_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01_earth2inst.nc');
            td_awac = rotate2(td_awac, 'earth');
            td_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo_earth2inst.nc');
            td_sig_ie = rotate2(td_sig_ie, 'earth');
            td_sig_o = rmfield(td_sig_o, 'orientmat');
            td_sig_o = rotate2(td_sig_o, 'earth');

            cd = read_netcdf('../../examples/data/dolfyn/control/RDI_test01_rotate_inst2earth.nc');
            cd_awac = read_netcdf('../../examples/data/dolfyn/control/AWAC_test01.nc');
            cd_sig_ie = read_netcdf('../../examples/data/dolfyn/control/Sig500_Echo.nc');
            cdwr2 = read_netcdf('../../examples/data/dolfyn/control/winriver02_rotate_ship2earth.nc');
            cd_sig = read_netcdf('../../examples/data/dolfyn/control/BenchFile01_rotate_inst2earth.nc');
            cd_sig_i = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU_rotate_inst2earth.nc');

            diff1 = Dolfyn_Test_Rotate.compare_structures(td, cd);
            diff2 = Dolfyn_Test_Rotate.compare_structures(tdwr2, cdwr2);
            diff3 = ...
                Dolfyn_Test_Rotate.compare_structures(td_awac,cd_awac);
            diff4 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig, cd_sig);
            diff5 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig_i, cd_sig_i);
            diff6 = ...
                Dolfyn_Test_Rotate.compare_structures(td_sig_ie,cd_sig_ie);

            tmp1 = isnan(cd_sig_ie.vel.data);
            tmp2 = isnan(td_sig_o.vel.data);
            tmp = tmp1|tmp2;

            dt1 = double(cd_sig_ie.vel.data);
            dt1(tmp) = 0.0;                        
            dt2 = double(td_sig_o.vel.data);
            dt2(tmp) = 0.0;  
            diff7 = abs(sum(abs(dt1 - dt2),...
                    1:numel(size(dt1)))/length(dt1));
            
            testCase.assertLessThan(diff1, 1e-5);
            testCase.assertLessThan(diff2, 1e-5);            
            testCase.assertLessThan(diff3, 1e-5);
            testCase.assertLessThan(diff4, 1e-5);
            testCase.assertLessThan(diff5, 1e-5);
            testCase.assertLessThan(diff6, 1e-5);
            testCase.assertLessThan(diff7, 1e-5);
        end

%         tr.dat_rdi = '../../examples/data/dolfyn/control/RDI_test01.nc';
%         tr.dat_sig = '../../examples/data/dolfyn/control/BenchFile01.nc';
%         tr.dat_sig_i = '../../examples/data/dolfyn/control/Sig1000_IMU.nc';
%         tr.dat_sig_ieb = '../../examples/data/dolfyn/control/VelEchoBT01.nc';
        
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
                fprintf('Field = %s | Diff = %f\n',field, diff)
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
                                ds_read.attrs.(field)));
                        else
                            diff = diff + sum(abs(ds_cntrl.attrs.(field)...
                                - ds_read.attrs.(field)), ...
                                1:numel(size(ds_cntrl.attrs.(field))))...
                                /length(ds_cntrl.attrs.(field));
                        end
                    end
                end
                fprintf('Field = %s | Diff = %f\n',field, diff)
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
                fprintf('Field = %s | Diff = %f\n',field, diff)
            end  
            %fprintf('Final Diff = %f\n',diff)
            format(oldFmt);
        end
    end
end

