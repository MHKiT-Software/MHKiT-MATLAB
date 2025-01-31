classdef Loads_TestLoads < matlab.unittest.TestCase

    methods (Test, TestTags = {'DebuggingActions'})

        function test_mhkit_import(testCase)
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

    end

    methods (Test)

        function test_bin_statistics(testCase)
            assumeFail(testCase, "Temporarily Disabled - Need to ask @hivanov about nan vs zero comparison")

            % create array containg wind speeds to use as bin edges
            bin_edges = 3:1:25;
            relative_file_name = '../../examples/data/loads/loads_data_dict.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string

            % Apply function to calculate means
            load_means_struct = [data.means.uWind_80m];
            load_means_table = struct2table(data.means);
            ta = struct2table([data.bin_means]);
            st = table2struct(ta,'ToScalar',true);
            ta_std = struct2table([data.bin_means_std]);
            st_std = table2struct(ta_std,'ToScalar',true);
            bin_against = load_means_struct;% load_means_table.uWind_80m;
            datast = bin_statistics(load_means_table,bin_against,bin_edges);
            assertEqual(testCase,st,datast.averages,'RelTol',0.001)
            assertEqual(testCase,st_std,datast.std,'RelTol',0.001)
        end

        function test_damage_equivalent_loads(testCase)
            relative_file_name = '../../examples/data/loads/loads_data_dict.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string

            fatigue_tower = 3804;
            fatigue_blade = 1388;

            loads_data_table = struct2table(data.loads);
            loads_data = table2struct(loads_data_table,'ToScalar',true);
            tower_load = loads_data.TB_ForeAft;
            blade_load = loads_data.BL1_FlapMom;
             DEL_tower = damage_equivalent_load(tower_load, 4,'bin_num',100,'data_length',600);
             DEL_blade = damage_equivalent_load(blade_load,10,'bin_num',100,'data_length',600);

             err_tower = abs((fatigue_tower-DEL_tower)/fatigue_tower);
             err_blade = abs((fatigue_blade-DEL_blade)/fatigue_tower);

             assertLessThan(testCase,err_tower,0.05)
             assertLessThan(testCase,err_blade,0.05)
        end

        function test_blade_moments(testCase)
            format long
            relative_file_name = '../../examples/data/loads/blade_cal.csv';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            blade_data = readmatrix(full_file_name);
            flap_offset = 9.19906E-05;
            edge_offset = -0.000310854;
            blade_matrix = [1034671.4,-126487.28,82507.959,1154090.7];

            flap_raw = blade_data(:,1);

            edge_raw = blade_data(:,2);

            [M_flap, M_edge] = blade_moments(blade_matrix,flap_offset,flap_raw,edge_offset,edge_raw);

            assertEqual(testCase,M_flap',blade_data(:,3),'RelTol',0.01);
            assertEqual(testCase,M_edge',blade_data(:,4),'RelTol',0.01);
        end

        function test_plot_statistics(testCase)
            relative_file_name = '../../examples/data/loads/loads_data_dict.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string

            means_data_table = struct2table(data.means);
            means_data = table2struct(means_data_table,'ToScalar',true);
            maxs_data_table = struct2table(data.maxs);
            maxs_data = table2struct(maxs_data_table,'ToScalar',true);
            mins_data_table = struct2table(data.mins);
            mins_data = table2struct(mins_data_table,'ToScalar',true);
            std_data_table = struct2table(data.std);
            std_data = table2struct(std_data_table,'ToScalar',true);

            % functionine path
            savepath = 'test_scatplotter.png';
            % Generate plot
            plot_statistics(means_data.uWind_80m,means_data.TB_ForeAft,maxs_data.TB_ForeAft,mins_data.TB_ForeAft,"y_stdev",std_data.TB_ForeAft,"xlabel",'Wind Speed [m/s]',"ylabel",'Tower Base Mom [kNm]',"savepath",savepath);

            assertTrue(testCase,isfile(savepath));
            delete(savepath);
        end

        function test_plot_bin_statistics(testCase)
            relative_file_name = '../../examples/data/loads/loads_data_dict.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string

            bin_means_data_table = struct2table(data.bin_means);
            bin_means_data = table2struct(bin_means_data_table,'ToScalar',true);
            bin_maxs_data_table = struct2table(data.bin_maxs);
            bin_maxs_data = table2struct(bin_maxs_data_table,'ToScalar',true);
            bin_mins_data_table = struct2table(data.bin_mins);
            bin_mins_data = table2struct(bin_mins_data_table,'ToScalar',true);

            bin_means_std_data_table = struct2table(data.bin_means_std);
            bin_means_std_data = table2struct(bin_means_std_data_table,'ToScalar',true);
            bin_maxs_std_data_table = struct2table(data.bin_maxs_std);
            bin_maxs_std_data = table2struct(bin_maxs_std_data_table,'ToScalar',true);
            bin_mins_std_data_table = struct2table(data.bin_mins_std);
            bin_mins_std_data = table2struct(bin_mins_std_data_table,'ToScalar',true);

            % functionine signal name, path, and bin centers
            savepath = 'test_binplotter.png';
%             bin_centers = np.arange(3.5,25.5,step=1)
            bin_centers = 3.5:1:24.5;
            signal_name = 'TB_ForeAft';

            % Specify inputs to be used in plotting
            bin_mean = bin_means_data.TB_ForeAft;
            bin_max  = bin_maxs_data.TB_ForeAft;
            bin_min  = bin_mins_data.TB_ForeAft;
            bin_mean_std = bin_means_std_data.TB_ForeAft;
            bin_max_std = bin_maxs_std_data.TB_ForeAft;
            bin_min_std = bin_mins_std_data.TB_ForeAft;

            % Generate plot
            plot_bin_statistics(bin_centers,bin_mean,bin_max,bin_min,bin_mean_std,bin_max_std,bin_min_std,"xlabel",'Wind Speed [m/s]',"ylabel",signal_name,"title",'Binned Stats',"savepath",savepath);

            assertTrue(testCase,isfile(savepath));
            delete(savepath);
        end

    end

end
