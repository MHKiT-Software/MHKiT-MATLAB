classdef Utils_TestGenUtils < matlab.unittest.TestCase

    methods (Test) 



        function test_get_statistics(testCase)
            relative_file_name = '../../examples/data/loads/loads_data_dict.json'; % filename in JSON extension
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string
            
            freq = 50; % Hz
            period = 600; % seconds
            vector_channels = {"WD_Nacelle","WD_NacelleMod"};
            
            % load in file
            loads_data_table = struct2table(data.loads);
            df = table2struct(loads_data_table,'ToScalar',true);
            
            df.Timestamp = datetime(df.Timestamp);
            df.time = df.Timestamp;
            % run function
            stats = get_statistics(df,freq,"period",period,"vector_channels",vector_channels);
            % check statistics
            assertEqual(testCase,stats.mean.uWind_80m,7.773,'AbsTol',0.01); % mean
            assertEqual(testCase,stats.max.uWind_80m,13.271,'AbsTol',0.01); % max
            assertEqual(testCase,stats.min.uWind_80m,3.221,'AbsTol',0.01); % min
            assertEqual(testCase,stats.std.uWind_80m,1.551,'AbsTol',0.01); % standard deviation3
            assertEqual(testCase,stats.std.WD_Nacelle,36.093,'AbsTol',0.01); % std vector averaging
            assertEqual(testCase,stats.mean.WD_Nacelle,178.1796,'AbsTol',0.01);% mean vector averaging
        end      

        function test_excel_to_datetime(testCase)
            % store excel timestamp
            excel_time = 42795.49212962963;
            % corresponding datetime
            time = datetime(2017,03,01,11,48,40);
            % test function
            answer = excel_to_datetime(excel_time);

            % check if answer is correct
            assertEqual(testCase,answer,time);
        end

        function test_magnitude_phase(testCase)
            % 2-d function
            magnitude = 9;
            y = sqrt(1/2*magnitude^2); x=y;
            phase = atan2(y,x);
            [mag, theta] = magnitude_phase({x; y});
            assert(all(magnitude == mag))
            assert(all(phase == theta))
            xx = [x,x]; yy = [y,y];
            [mag, theta] = magnitude_phase({xx; yy});
            assert(all(magnitude == mag))
            assert(all(phase == theta))
            % 3-d function
            magnitude = 9;
            y = sqrt(1/3*magnitude^2); x=y; z=y;
            phase1 = atan2(y,x);
            phase2 = atan2(sqrt(x.^2 + y.^2),z);
            [mag, theta, phi] = magnitude_phase({x; y; z});
            assert(all(magnitude == mag))
            assert(all(phase1 == theta))
            assert(all(phase2 == phi))
            xx = [x,x]; yy = [y,y]; zz = [z,z];
            [mag, theta, phi] = magnitude_phase({xx; yy; zz});
            assert(all(magnitude == mag))
            assert(all(phase1 == theta))
            assert(all(phase2 == phi))
        end
        function test_read_nc_file_group(testCase)
            % MATLAB version should >= 2021b
            
            %1. Check LongName with group path
            fnm = 'QA4ECV_L2_NO2_OMI_20180301T052400_o72477_fitB_v1.nc';
            res = read_nc_file(strcat('example_ncfiles/',fnm));
            val1 = res.groups.PRODUCT.groups.SUPPORT_DATA.groups.INPUT_DATA.LongName;
            val2 = '/PRODUCT/SUPPORT_DATA/INPUT_DATA';
            assertEqual(testCase,val1,val2);
            %2. Check Group Attributes
            finfo = ncinfo(strcat('example_ncfiles/',fnm));
            val1 = res.groups.METADATA.groups.ALGORITHM_SETTINGS.groups.SLANT_COLUMN_RETRIEVAL.Attributes;
            val2 = finfo.Groups(2).Groups(1).Groups(1).Attributes;
            assertEqual(testCase,val1,val2);
            %3. Check Variables: file with groups
            % '/PRODUCT/SUPPORT_DATA/DETAILED_RESULTS'
            ginfo = finfo.Groups(1).Groups(1).Groups(2);
            vnms = {ginfo.Variables.Name};
            sz = size(ginfo.Variables);
            % 3.1 check Dims
            idx = randi([1,sz(2)],1); 
            vname = check_name(vnms{idx});
            val1 = res.groups.PRODUCT.groups.(['SUPPORT_' ...
                'DATA']).groups.DETAILED_RESULTS.Variables.(vname).Dims;
            val2 = {ginfo.Variables(idx).Dimensions.Name};
            val3 = size(res.groups.PRODUCT.groups.(['SUPPORT_' ...
                'DATA']).groups.DETAILED_RESULTS.Variables.(vname).Data);
            val4 = size(ncread(strcat('example_ncfiles/',fnm),...
                strcat('PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/',vnms{idx})));
            assertEqual(testCase,val1,val2);
            assertEqual(testCase,val3,val4);
            % 3.2 check Data
            idx = randi([1,sz(2)],1); 
            vname = check_name(vnms{idx});
            val1 = res.groups.PRODUCT.groups.(['SUPPORT_' ...
                'DATA']).groups.DETAILED_RESULTS.Variables.(vname).Data;
            val2 = ncread(strcat('example_ncfiles/',fnm),...
                strcat('PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/',vnms{idx}));
            
            testCase.verifyTrue(isequaln(val1,val2),vname);
            % 3.3 check Attributes
            idx = randi([1,sz(2)],1); 
            vname = check_name(vnms{idx});
            in_names = {ginfo.Variables(idx).Attributes.Name};
            in_vals = {ginfo.Variables(idx).Attributes.Value};
            xtemp = res.groups.PRODUCT.groups.(['SUPPORT_' ...
                'DATA']).groups.DETAILED_RESULTS.Variables;
            for iattr = 1:numel(in_names)
                if strcmp(in_names{iattr},'_FillValue')
                    out_val = xtemp.(vname).FillValue;
                else
                    out_val = xtemp.(vname).Attrs.(in_names{iattr});
                end
                testCase.verifyTrue(isequaln(in_vals{iattr},out_val));
            end
        end
        
        function test_read_nc_file_nogroup(testCase)
            % MATLAB version should >= 2021b
            fnms = {dir('example_ncfiles/').name};
            % file without groups: check variables
            for ifnm = 3:numel(fnms)
                fnm = fnms{ifnm};
                if strcmp(fnm,'QA4ECV_L2_NO2_OMI_20180301T052400_o72477_fitB_v1.nc')
                    continue
                end
                fprintf("Checking File: %s \n",fnm);
                res = read_nc_file(strcat('example_ncfiles/',fnm));
                ginfo = ncinfo(strcat('example_ncfiles/',fnm));
                vnms = {ginfo.Variables.Name};
                sz = size(ginfo.Variables);
                % 1 check Dims
                idx = randi([1,sz(2)],1); 
                count = 0;
                while (isempty(ginfo.Variables(idx).Dimensions)&&count<10)
                    idx = randi([1,sz(2)],1); 
                    count = count + 1;
                end
                if ~isempty(ginfo.Variables(idx).Dimensions)
                    vname = check_name(vnms{idx});
                    val1 = res.Variables.(vname).Dims;
                    val2 = {ginfo.Variables(idx).Dimensions.Name};
                    val3 = size(res.Variables.(vname).Data);
                    val4 = size(ncread(strcat('example_ncfiles/',fnm),vnms{idx}));
                    assertEqual(testCase,val1,val2);
                    assertEqual(testCase,val3,val4);
                end
                % 2 check Data
                idx = randi([1,sz(2)],1); 
                count = 0;
                while (isempty(ginfo.Variables(idx).Dimensions)&&count<10)
                    idx = randi([1,sz(2)],1); 
                    count = count + 1;
                end
                vname = check_name(vnms{idx});
                val1 = res.Variables.(vname).Data;
                val2 = ncread(strcat('example_ncfiles/',fnm),vnms{idx});
                testCase.verifyTrue(isequaln(val1,val2),vname);
                % 3 check Attributes
                idx = randi([1,sz(2)],1); 
                count = 0;
                while (isempty(ginfo.Variables(idx).Attributes)&&count<10)
                    idx = randi([1,sz(2)],1); 
                    count = count + 1;
                end
                vname = check_name(vnms{idx});
                if ~isempty(ginfo.Variables(idx).Attributes)
                    in_names = {ginfo.Variables(idx).Attributes.Name};
                    in_vals = {ginfo.Variables(idx).Attributes.Value};
                    for iattr = 1:numel(in_names)
                        if strcmp(in_names{iattr},'_FillValue')
                            out_val = res.Variables.(vname).FillValue;
                        else
                            out_val = res.Variables.(vname).Attrs.(in_names{iattr});
                        end
                        testCase.verifyTrue(isequaln(in_vals{iattr},out_val));
                    end
                end
            
            end
        end
        
        function test_read_nc_file_var(testCase)
            % MATLAB version should >= 2021b
            fnms = {dir('example_ncfiles/').name};
            %1. test on file with group:
            fnm = 'QA4ECV_L2_NO2_OMI_20180301T052400_o72477_fitB_v1.nc';
            varlst = {'PRODUCT/amf_total','PRODUCT/amf_trop',...
                'PRODUCT/latitude','PRODUCT/averaging_kernel',...
                'PRODUCT/tropospheric_no2_vertical_column',...
                'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_zenith_angle',...
                'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds',...
                'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/radiance_calibration_stretch'};
            res = read_nc_file_var(strcat('example_ncfiles/',fnm),...
                varlst,0);
            res1 = read_nc_file_var(strcat('example_ncfiles/',fnm),...
                varlst,1);           
            idx = randi([1,length(varlst)],1);
            var2check = varlst{idx}; nstr = split(var2check,'/');
            vname = check_name(nstr{end});
            % 1.1 check Data Field:
            val1 = res.(vname).Data;
            val1_1 = res1(idx).Data;
            val2 = ncread(strcat('example_ncfiles/',fnm),var2check);
            %1.1.1 test opt=0: output as struct
            testCase.verifyTrue(isequaln(val1,val2),...
                strcat('opt=0,',var2check));
            %1.1.1 test opt=1: output as struct array
            testCase.verifyTrue(isequaln(val1_1,val2),...
                strcat('opt=1,',var2check));
            % 1.2 check Dims Names:
            vinfo = ncinfo(strcat('example_ncfiles/',fnm),var2check);
            val1 = res.(vname).Dims;
            val1_1 = res1(idx).Dims;
            val2 = {vinfo.Dimensions.Name};
            testCase.verifyTrue(isequaln(val1,val2),...
                strcat('opt=0,',var2check));
            testCase.verifyTrue(isequaln(val1_1,val2),...
                strcat('opt=1,',var2check));
            %1.3 check Attrs & FillValue:
            val1 = res.(vname).FillValue;
            val1_1 = res1(idx).FillValue;
            val2 = vinfo.FillValue;
            testCase.verifyTrue(isequaln(val1,val2),...
                strcat('opt=0,',var2check));
            testCase.verifyTrue(isequaln(val1_1,val2),...
                strcat('opt=1,',var2check));
            if ~isempty(vinfo.Attributes)
                in_names = {vinfo.Attributes.Name};
                in_vals = {vinfo.Attributes.Value};
                for iattr = 1:numel(in_names)
                    if strcmp(in_names{iattr},'_FillValue')
                        out_val = res.(vname).FillValue;
                        out_val1 = res1(idx).FillValue;
                    else
                        out_val = res.(vname).Attrs.(in_names{iattr});
                        out_val1 = res1(idx).Attrs(iattr).Value;
                    end
                    testCase.verifyTrue(isequaln(in_vals{iattr},...
                        out_val),strcat('opt=0,',var2check,': ',in_names{iattr}));
                    testCase.verifyTrue(isequaln(in_vals{iattr},...
                        out_val1),strcat('opt=1,',var2check,': ',in_names{iattr}));
                end
            end
            
            %2. test on files without group:
            for ifnm = 3:numel(fnms)
                fnm = fnms{ifnm};
                if strcmp(fnm,'QA4ECV_L2_NO2_OMI_20180301T052400_o72477_fitB_v1.nc')
                    continue
                end
                fprintf("Checking File: %s \n",fnm);
                ginfo = ncinfo(strcat('example_ncfiles/',fnm));                
                vnms = {ginfo.Variables.Name};
                res = read_nc_file_var(strcat('example_ncfiles/',fnm),...
                    vnms,0);
                res1 = read_nc_file_var(strcat('example_ncfiles/',fnm),...
                    vnms,1);
                sz = length(vnms);
                % 2.1 check Data Field: 
                idx = randi([1,sz],1);%max([fix(sz*0.1),1]));
                var2check = vnms{idx};
                vname = check_name(var2check);
                val1 = res.(vname).Data;
                val1_1 = res1(idx).Data;
                val2 = ncread(strcat('example_ncfiles/',fnm),var2check);
                %2.1.1 test opt=0: output as struct
                testCase.verifyTrue(isequaln(val1,val2),...
                    strcat('opt=0,',var2check));
                %2.1.1 test opt=1: output as struct array
                testCase.verifyTrue(isequaln(val1_1,val2),...
                    strcat('opt=1,',var2check));
                %testCase.verifyTrue(isequaln(val1,val2),var2check);
                % 2.2 check Dims Names:
                vinfo = ncinfo(strcat('example_ncfiles/',fnm),var2check);
                if ~isempty(vinfo.Dimensions)
                    val1 = res.(vname).Dims;
                    val1_1 = res1(idx).Dims;
                    val2 = {vinfo.Dimensions.Name};
                    testCase.verifyTrue(isequaln(val1,val2),...
                        strcat('opt=0,',var2check));
                    testCase.verifyTrue(isequaln(val1_1,val2),...
                        strcat('opt=1,',var2check));
                    %testCase.verifyTrue(isequaln(val1,val2),var2check);
                else
                    testCase.verifyTrue(isempty(res.(vname).Dims),...
                        strcat('opt=0,',var2check));
                    testCase.verifyTrue(isempty(res1(idx).Dims),...
                        strcat('opt=1,',var2check));
                end
                % 2.3 check Attrs & FillValue:
                val1 = res.(vname).FillValue;
                val1_1 = res1(idx).FillValue;
                val2 = vinfo.FillValue;
                testCase.verifyTrue(isequaln(val1,val2),...
                    strcat('opt=0,',var2check));
                testCase.verifyTrue(isequaln(val1_1,val2),...
                    strcat('opt=1,',var2check));
                %testCase.verifyTrue(isequaln(val1,val2),var2check);
                if ~isempty(vinfo.Attributes)
                    in_names = {vinfo.Attributes.Name};
                    in_vals = {vinfo.Attributes.Value};
                    for iattr = 1:numel(in_names)
                        if strcmp(in_names{iattr},'_FillValue')
                            out_val = res.(vname).FillValue;
                            out_val_1 = res1(idx).FillValue;
                        else
                            out_val = res.(vname).Attrs.(in_names{iattr});
                            out_val_1 = res1(idx).Attrs(iattr).Value;
                        end
                        testCase.verifyTrue(isequaln(in_vals{iattr},...
                            out_val),strcat('opt=0,',var2check,': ',in_names{iattr}));
                        testCase.verifyTrue(isequaln(in_vals{iattr},...
                            out_val_1),strcat('opt=1,',var2check,': ',in_names{iattr}));
                    end
                end
            end
            
            
        end
    end
end  




