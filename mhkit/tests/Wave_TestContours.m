classdef Wave_TestContours < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function test_samples_contour(testCase)
            file_loc = "../../examples/data/wave/WDRT_caluculated_countours.json";
            str = fileread(file_loc); % dedicated for reading files as text 
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string 
            te_samples = [10, 15, 20];
            hs_samples_0 = [8.56637939, 9.27612515, 8.70427774];
            hs_contour = data.gaussian_x1;
            te_contour = data.gaussian_x2;
            hs_samples = samples_contour(te_samples, te_contour, hs_contour);
            
            assertEqual(testCase, hs_samples, hs_samples_0, 'AbsTol', 0.0005)
        end

        function test_samples_full_seastate(testCase)
            hs_0 = [5.91760129, 4.55185088, 1.41144991, 12.64443154, 7.89753791, 0.93890797];
            te_0 = [14.24199604, 8.25383556, 6.03901866, 16.9836369, 9.51967777, 3.46969355];
            w_0 = [2.18127398e-01,2.18127398e-01,2.18127398e-01,2.45437862e-07,2.45437862e-07,2.45437862e-07];

            file_loc = "../../examples/data/wave/Hm0_Te_46022.json";
            str = fileread(file_loc); % dedicated for reading files as text 
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Hm0 = cell2mat(struct2cell(data.Hm0));
            qc = find(Hm0 < 20);
            Te = cell2mat(struct2cell(data.Te));
            dt_ss = 3600;
            points_per_interval = 3;
            return_periods = [50, 100];
            py.numpy.random.seed(int8(0));
            [hs, te, w] = samples_full_seastate(Hm0(qc), Te(qc), points_per_interval, return_periods, dt_ss, "PCA", 250);
            assertEqual(testCase, hs, hs_0, 'AbsTol',0.0005)
            assertEqual(testCase, te, te_0, 'AbsTol',0.0005)
            assertEqual(testCase, w, w_0, 'AbsTol',0.0005)
        end

        function test_environmental_contour(testCase)

            % assumeFail(testCase, "Not compatible with latest MHKIT-Python")

            relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);

            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Te_table = struct2table(valdata1.Te,'AsArray',true);
            Te = table2array(Te_table);
            Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
            Hm0 = table2array(Hm0_table);

            filter = Hm0 < 20;
            Hm0 = Hm0(filter);
            Te = Te(filter);
            [row, col] = find(~isnan(Te));
            Hm0 = Hm0(col);
            Te = Te(col);
            [row, col] = find(~isnan(Hm0));
            Hm0 = Hm0(col);
            Te = Te(col);

            time_str = Hm0_table.Properties.VariableNames;

            time1 = str2num(erase(time_str{1},'x'));
            time2 = str2num(erase(time_str{2},'x'));

            dt = (time2-time1)/1000.;
            time_R = 100;

            contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');

            relative_file_name= '../../examples/data/wave/Hm0_Te_contours_46022.csv';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            expected_contours = readmatrix(full_file_name);

            Hm0_expected = expected_contours(:,1);
            Te_expected = expected_contours(:,2);

            assertEqual(testCase,contour.contour1,Hm0_expected','RelTol',0.01);
            assertEqual(testCase,contour.contour2,Te_expected','RelTol',0.01);


        end

        function test_plot_environmental_contour(testCase)

            %assumeFail(testCase, "Not compatible with latest MHKIT-Python")

            relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);

            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Te_table = struct2table(valdata1.Te,'AsArray',true);
            Te = table2array(Te_table);
            Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
            Hm0 = table2array(Hm0_table);

            filter = Hm0 < 20;
            Hm0 = Hm0(filter);
            Te = Te(filter);
            [row, col] = find(~isnan(Te));
            Hm0 = Hm0(col);
            Te = Te(col);
            [row, col] = find(~isnan(Hm0));
            Hm0 = Hm0(col);
            Te = Te(col);

            time_str = Hm0_table.Properties.VariableNames;

            time1 = str2num(erase(time_str{1},'x'));
            time2 = str2num(erase(time_str{2},'x'));

            dt = (time2-time1)/1000.;
            time_R = 100;

            contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');

            filename = 'wave_plot_env_contour.png';
            if isfile(filename)
                delete(filename);
            end


            plot_environmental_contours(Te, Hm0,contour.contour2,contour.contour1,"savepath",filename...
                ,"x_label",...
                'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
                "contour_label",'100 Year Contour');
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        % function test_plot_environmental_contour_multiyear(testCase)
        %
        %     assumeFail(testCase, "Not compatible with latest MHKIT-Python")
        %     % not sure about why this test exists...return period has to be float or
        %     % int and cannot be a list...
        %     relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
        %     full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
        %
        %     fid = fopen(full_file_name); % Opening the file
        %     raw = fread(fid,inf); % Reading the contents
        %     str = char(raw'); % Transformation
        %     fclose(fid); % Closing the file
        %     valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
        %     Te_table = struct2table(valdata1.Te,'AsArray',true);
        %     Te = table2array(Te_table);
        %     Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
        %     Hm0 = table2array(Hm0_table);
        %
        %     filter = Hm0 < 20;
        %     Hm0 = Hm0(filter);
        %     Te = Te(filter);
        %     [row, col] = find(~isnan(Te));
        %     Hm0 = Hm0(col);
        %     Te = Te(col);
        %     [row, col] = find(~isnan(Hm0));
        %     Hm0 = Hm0(col);
        %     Te = Te(col);
        %
        %     time_str = Hm0_table.Properties.VariableNames;
        %
        %     time1 = str2num(erase(time_str{1},'x'));
        %     time2 = str2num(erase(time_str{2},'x'));
        %
        %     dt = (time2-time1)/1000.;
        %     time_R = [100, 120, 130];
        %
        %     contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');
        %
        %     filename = 'wave_plot_env_contour_multiyear.png';
        %     if isfile(filename)
        %         delete(filename);
        %     end
        %
        %
        %     plot_environmental_contours(Te, Hm0,contour.contour2,contour.contour1,"savepath",filename...
        %         ,"x_label",...
        %         'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
        %         "contour_label",{'100 Year Contour','120 Year Contour','130 Year Contour'});
        %     assertTrue(testCase,isfile(filename));
        %     delete(filename);
        % end

    end
    
end