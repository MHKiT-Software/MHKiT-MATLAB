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
                % Only check numeric variables. For string and char variables
                % (e.g. inst, earth, dir) read_nc_file_var reports a NaN
                % FillValue while ncinfo reports "", so they always fail.
                numeric_idx = find(~ismember({ginfo.Variables.Datatype}, {'string', 'char'}));
                % 2.1 check Data Field:
                idx = numeric_idx(randi(numel(numeric_idx)));
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

        function test_mhkit_standardize_user_input_to_column_vectors(testCase)
            % Row vector is transposed and flagged
            row = [1 2 3 4];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(row, 'function_name', 'test_fn');
            assertEqual(testCase, out, row(:));
            assertTrue(testCase, was_row);

            % Column vector passes through unchanged, not flagged
            col = [1; 2; 3; 4];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(col, 'function_name', 'test_fn');
            assertEqual(testCase, out, col);
            assertFalse(testCase, was_row);

            % Scalar passes through unchanged, not flagged
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(5, 'function_name', 'test_fn');
            assertEqual(testCase, out, 5);
            assertFalse(testCase, was_row);

            % Matrix (multiple column-oriented vectors) passes through unchanged
            mat = [1 2; 3 4; 5 6];
            [out, was_row] = mhkit_standardize_user_input_to_column_vectors(mat, 'function_name', 'test_fn');
            assertEqual(testCase, out, mat);
            assertFalse(testCase, was_row);

            % Empty input errors with the calling function's name in the identifier
            testCase.verifyError(@() mhkit_standardize_user_input_to_column_vectors([], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');

            % 3-D input errors with the calling function's name in the identifier
            testCase.verifyError(@() mhkit_standardize_user_input_to_column_vectors(ones(2,2,2), 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');
        end

        function test_mhkit_restore_column_vectors_to_user_input(testCase)
            col = [1; 2; 3; 4];

            % was_row = true restores to a row vector
            out = mhkit_restore_column_vectors_to_user_input(col, true);
            assertEqual(testCase, out, col.');

            % was_row = false leaves the column vector unchanged
            out = mhkit_restore_column_vectors_to_user_input(col, false);
            assertEqual(testCase, out, col);

            % Round trip through standardize + restore recovers the original orientation
            row = [1 2 3 4];
            [standardized, was_row] = mhkit_standardize_user_input_to_column_vectors(row, 'function_name', 'test_fn');
            restored = mhkit_restore_column_vectors_to_user_input(standardized, was_row);
            assertEqual(testCase, restored, row);
        end

        function test_mhkit_verify_is_column_vector(testCase)
            % Column vector and scalar pass silently, numeric or datetime/duration
            mhkit_verify_is_column_vector([1;2;3], 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(5, 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(datetime(2026,1,1) + hours(0:1)', 'function_name', 'test_fn');
            mhkit_verify_is_column_vector(hours(0:1)', 'function_name', 'test_fn');

            % Row vector and matrix both raise a function-scoped error
            testCase.verifyError(@() mhkit_verify_is_column_vector([1 2 3], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidOutput');
            testCase.verifyError(@() mhkit_verify_is_column_vector([1 2; 3 4], 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidOutput');

            % Unsupported type raises a function-scoped error
            testCase.verifyError(@() mhkit_verify_is_column_vector("not numeric", 'function_name', 'my_fn'), ...
                'MHKiT:my_fn:InvalidInput');
        end

        function test_mhkit_frequency_to_column_names(testCase)
            frequency = [0.1; 0.2; 0.3; 0.4];
            names = mhkit_frequency_to_column_names(frequency);
            assertEqual(testCase, names, ["f_0_1000Hz"; "f_0_2000Hz"; "f_0_3000Hz"; "f_0_4000Hz"]);

            % Every generated name must be a valid MATLAB identifier
            for i = 1:numel(names)
                assertTrue(testCase, isvarname(names(i)));
            end

            % Custom decimal precision
            names2 = mhkit_frequency_to_column_names([0.1; 0.2], 'decimals', 2);
            assertEqual(testCase, names2, ["f_0_10Hz"; "f_0_20Hz"]);

            % Frequencies that collide once formatted must error, not
            % silently produce duplicate table variable names
            testCase.verifyError(@() mhkit_frequency_to_column_names([0.101; 0.102], 'decimals', 2), ...
                'MHKiT:mhkit_frequency_to_column_names:DuplicateNames');
        end

    end

end

