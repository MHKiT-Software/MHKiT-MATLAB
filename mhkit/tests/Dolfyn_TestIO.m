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
            diff = 0.0;
            exclude = {'coords', 'attrs', 'time', 'complex_vars'};
            % Check coords first
            fields = fieldnames(ds_read.coords);
            for qq = 1:numel(fields)
                field = fields{qq}
                if ~any(strcmp(exclude,field))
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
                            diff = (diff + ...
                                sum(abs(double(ds_cntrl.coords.(field)) -...
                                double(ds_read.coords.(field)'))))/2.0;
                        else
                            diff = (diff + ...
                                sum(abs(double(ds_cntrl.coords.(field)) -...
                                double(ds_read.coords.(field)))))/2.0;
                        end
                    end
                end
                fprintf('Diff = %f\n',diff)
            end
    
            % Check Attributes
            fields = fieldnames(ds_cntrl.attrs);
            for qq = 1:numel(fields)
                field = fields{qq};            
                if ~any(strcmp(exclude,field))
                    if iscell(ds_cntrl.attrs.(field))                
                        for kk = 1:numel(ds_cntrl.attrs.(field))
                            chk_nm = ds_cntrl.attrs.(field){kk};
                            diff = diff + ...
                                double(~any(strcmpi(chk_nm,...
                                ds_read.attrs.(field){:})));
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
                    elseif islogical(ds_cntrl.attrs.(field))
                        debug = 1;
                    elseif isnumeric(ds_cntrl.attrs.(field))
                        diff = diff + ... 
                            sum(double(ds_cntrl.attrs.(field) ~=...
                            ds_read.attrs.(field)'));
                    end
                end
                fprintf('Diff = %f\n',diff)
            end
    
            % Now check the remaining fields
            fields = fieldnames(ds_cntrl);
            for qq = 1:numel(fields)
                field = fields{qq}
                if strcmpi(field,'vel_bt')
                    degug = 1;
                end
                if ~any(strcmp(exclude,field))
                    cls = class(ds_cntrl.(field));
                    if strcmp(cls,'struct')
                        % Loop through structure
                        % Data
                        tmp1 = isnan(ds_cntrl.(field).data);
                        dt1 = ds_cntrl.(field).data;
                        dt1(tmp1) = 0.0;
                        tmp2 = isnan(ds_read.(field).data);
                        dt2 = ds_read.(field).data;
                        dt2(tmp2) = 0.0;
                        if strcmpi(field,'orientmat')
                            % The difference between significant figures in pi,
                            % cos, sin, etc between numpy and matlab causes
                            % an inherent difference in the answers for
                            % orientmat
                            diff = (diff + ...
                                    sum(abs(single(round(dt1,4))...
                                    - single(round(dt1,4))),...
                                    1:numel(size(dt1))))/2.0;
                        else
                            diff = (diff + ...
                                    sum(abs(single(dt1) - single(dt2)),...
                                    1:numel(size(dt1))))/2.0;
                        end
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
                fprintf('Diff = %f\n',diff)
            end  
            fprintf('Final Diff = %f\n',diff)
        end
    end
end

