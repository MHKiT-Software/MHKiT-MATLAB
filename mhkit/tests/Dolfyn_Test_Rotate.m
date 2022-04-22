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

        function test_rotate_inst2earth(testCase)
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
                                ds_read.attrs.(field)'));
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
                fprintf('Field = %s | Diff = %f\n',field, diff)
            end  
            %fprintf('Final Diff = %f\n',diff)
            format(oldFmt);
        end
    end
end

