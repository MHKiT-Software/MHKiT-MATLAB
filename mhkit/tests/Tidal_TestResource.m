classdef Tidal_TestResource < matlab.unittest.TestCase

    methods (Test)

        function test_exceedance_probability(testCase)
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8], 'time', [0 1 2 3 4 5 6 7 8]);
            f = exceedance_probability(Q);
            assertEqual(testCase,min(f.F), 10);
            assertEqual(testCase,max(f.F), 90);
        end

        function test_principal_flow_directions(testCase)
            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;

            width_direction = 10;
            [direction1, direction2] = principal_flow_directions(data.d, width_direction);
            assertEqual(testCase,direction1,172.0);
            assertEqual(testCase,round(direction2,1),round(352.3,1));
        end

        function test_plot_current_timeseries(testCase)
            filename = 'tidal_plot_current_timeseries.png';
            if isfile(filename)
                delete(filename);
            end

            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;
            width_direction = 10;
            [direction1, direction2] = principal_flow_directions(data.d, width_direction);

            plot_current_timeseries(data,direction1,"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_plot_joint_probability_distribution(testCase)
            filename = 'tidal_plot_joint_probability_distribution.png';
            if isfile(filename)
                delete(filename)
            end

            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;
            width_direction = 10;

            data = rmfield(data,'id');
            data = rmfield(data,'name');
            data = rmfield(data,'lat');
            data = rmfield(data,'lon');
            data = rmfield(data,'b');

            plot_joint_probability_distribution(data,width_direction,0.1,"savepath",filename);

            assertTrue(testCase,isfile(filename))
            delete(filename);
        end

        function test_plot_rose(testCase)
            filename = 'tidal_plot_rose.png';
            if isfile(filename)
                delete(filename);
            end

            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;
            width_direction = 10;

            plot_rose(data,width_direction,1,"savepath",filename);

            assertTrue(testCase,isfile(filename))
            delete(filename);
        end

        function test_plot_phase_probability(testCase)
            filename = 'tidal_plot_phase_probability.png';
            if isfile(filename)
                delete(filename);
            end

            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;
            width_direction = 1;
            [direction1, direction2] = ...
                principal_flow_directions(data.d,width_direction);
            flood = direction1 ;
            ebb = direction2 ;

            plot_tidal_phase_probability(data,flood,ebb,"savepath",filename);

            assertTrue(testCase,isfile(filename))
            delete(filename);
        end

        function test_plot_phase_exceedance(testCase)

            assumeFail(testCase, "TODO: Fix - Error using area X must be same length as Y.");

            filename = 'tidal_plot_phase_exceedance.png';
            if isfile(filename)
                delete(filename);
            end

            relative_file_name = '../../examples/data/tidal/s08010.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_noaa_json(full_file_name);
            data.s = data.s/100;
            width_direction = 1;
            [direction1, direction2] = ...
                principal_flow_directions(data.d,width_direction);
            flood = direction1 ;
            ebb = direction2 ;

            plot_tidal_phase_exceedance(data,flood,ebb,"savepath",filename);

            assertTrue(testCase,isfile(filename))
            delete(filename);
        end

    end

end

