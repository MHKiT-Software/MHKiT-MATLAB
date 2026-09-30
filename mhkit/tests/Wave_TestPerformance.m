classdef Wave_TestPerformance < matlab.unittest.TestCase

    methods (Test)

        function test_capture_width(testCase)
            a = 40;
            b = 200;
            Obj.P = (b-a).*rand(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*rand(1,100000) + a;
            %Obj.J = normrnd(300, 10, [1,100000]);

            CW = capture_width(Obj.P, Obj.J);
            CW_stats = mean(CW);
            assertEqual(testCase,CW_stats, 1.4, 'RelTol',0.1);
        end

        function test_capture_width_matrix(testCase)
            seednum = 123;
            rng(seednum);
            a = 0.8;
            b = 4.5;
            Obj.Te = (b-a).*randn(1,100000) + a;
            %Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;
            %Obj.J = normrnd(300, 10, [1,100000]);
            sigma = 4;
            Obj.Hm0 = abs(sigma*randn(1,100000)+1i*sigma*randn(1,100000));
            %Obj.Hm0 = raylrnd(4, [1,100000]);
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;

            CW = capture_width(Obj.P, Obj.J);
            CWM = capture_width_matrix(Obj.Hm0, Obj.Te, CW, 'std', Obj.Hm0_bins, Obj.Te_bins);

            assertEqual(testCase,size(CWM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(CWM.values))), 43);
        end

        function test_wave_energy_flux_matrix(testCase)
            seednum = 123;
            rng(seednum);
            a = 0.8;
            b = 4.5;
            Obj.Te = (b-a).*randn(1,100000) + a;
            %Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;
            sigma = 4;
            Obj.Hm0 = abs(sigma*randn(1,100000)+1i*sigma*randn(1,100000));
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;

            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            assertEqual(testCase,size(JM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(JM.values))), 43);
        end

        function test_power_matrix(testCase)
            seednum = 123;
            rng(seednum);
            a = 0.8;
            b = 4.5;
            Obj.Te = (b-a).*randn(1,100000) + a;
            %Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;
            sigma = 4;
            Obj.Hm0 = abs(sigma*randn(1,100000)+1i*sigma*randn(1,100000));
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;

            CW = capture_width(Obj.P, Obj.J);
            CWM = capture_width_matrix(Obj.Hm0, Obj.Te,CW, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            PM = power_matrix(CWM, JM);
            assertEqual(testCase,size(PM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(PM.values))), 43);
        end

        function test_mean_annual_energy_production(testCase)
            rng(123);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;

            CW = capture_width(Obj.P, Obj.J);
            maep = mean_annual_energy_production_timeseries(CW, Obj.J);
            % CW .* J is P, so MAEP is the mean power over 8766 hours
            expected = 8766 * mean(Obj.P);
            assertEqual(testCase, maep, expected, 'RelTol', 1e-10);
        end

        function test_plot_matrix(testCase)
            filename = 'wave_plot_matrix.png';
            if isfile(filename)
                delete(filename);
            end

            seednum = 123;
            rng(seednum);
            a = 0.8;
            b = 4.5;
            Obj.Te = (b-a).*randn(1,100000) + a;
            %Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;
            sigma = 4;
            Obj.Hm0 = abs(sigma*randn(1,100000)+1i*sigma*randn(1,100000));
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            M = wave_energy_flux_matrix(Obj.Hm0,Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);

            plot_matrix(M,'Wave Energy Flux Matrix',"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_power_performance_workflow(testCase)
            filename = 'Capture Width Matrix mean.png';
            if isfile(filename)
                delete(filename);
            end

            relative_file_name = './../../examples/data/wave/data.txt';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            S1 = read_NDBC_file(full_file_name);
            h = 60;
            seednum = 123;
            rng(seednum);
            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,743) + a;

            [x, y] = power_performance_workflow(S1,h,Obj.P,"mean","savepath",'./');

            assertTrue(testCase,isfile(filename));
            delete(filename);
            % Only the requested statistics plus mean and probability (for MAEP)
            assertTrue(testCase,isfield(x,"mean"));
            assertTrue(testCase,isfield(x,"probability"));
            assertFalse(testCase,isfield(x,"std"));
            assertEqual(testCase,y,401239.4822345051, 'RelTol',0.00001);

        end

    end

end
