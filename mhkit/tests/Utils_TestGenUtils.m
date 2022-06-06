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
    end
end  




        