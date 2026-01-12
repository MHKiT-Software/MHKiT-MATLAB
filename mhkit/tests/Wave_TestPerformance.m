classdef Wave_TestPerformance < matlab.unittest.TestCase

    methods (Static)
        function bins = make_bins(bin_start, bin_end, bin_interval)
            % Generate bins for Python's capture_length_matrix.
            % Python expects "bin centers" but internally converts them to
            % edges. To maintain compatibility, we pass the edge values
            % directly as "centers" (excluding the final edge).
            bins.edges = bin_start:bin_interval:bin_end;
            bins.centers = bins.edges(1:end-1);  % Pass to Python as "centers"
        end
    end

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

            L = capture_length(Obj.P, Obj.J);
            L_stats = mean(L);
            assertEqual(testCase,L_stats, 1.4, 'RelTol',0.1);
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
            Hm0_bins = Wave_TestPerformance.make_bins(0, 19, 0.5);
            Te_bins = Wave_TestPerformance.make_bins(0, 9, 1);

            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te, L, 'std', Hm0_bins.centers, Te_bins.centers);

            assertEqual(testCase,size(LM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(LM.values))), 45, 'RelTol',0.1);
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
            Hm0_bins = Wave_TestPerformance.make_bins(0, 19, 0.5);
            Te_bins = Wave_TestPerformance.make_bins(0, 9, 1);

            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Hm0_bins.centers, Te_bins.centers);
            assertEqual(testCase,size(JM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(JM.values))), 45, 'RelTol',0.1);
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
            Hm0_bins = Wave_TestPerformance.make_bins(0, 19, 0.5);
            Te_bins = Wave_TestPerformance.make_bins(0, 9, 1);

            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te,L, 'mean', Hm0_bins.centers, Te_bins.centers);
            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Hm0_bins.centers, Te_bins.centers);
            PM = power_matrix(LM, JM);
            assertEqual(testCase,size(PM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(PM.values))), 45, 'RelTol',0.1);
        end

        function test_mean_annual_energy_production(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.AbsoluteTolerance

            a = 40;
            b = 200;
            Obj.P = (b-a).*randn(1,100000) + a;
            %Obj.P = normrnd(200, 40, [1,100000]);
            a = 10;
            b = 300;
            Obj.J = (b-a).*randn(1,100000) + a;

            L = capture_length(Obj.P, Obj.J);
            maep = mean_annual_energy_production_timeseries(L, Obj.J);
            testCase.verifyThat(1754020.077,IsEqualTo(maep,'Within',AbsoluteTolerance(2e+06)))
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
            Hm0_bins = Wave_TestPerformance.make_bins(0, 19, 0.5);
            Te_bins = Wave_TestPerformance.make_bins(0, 9, 1);
            M = wave_energy_flux_matrix(Obj.Hm0,Obj.Te,Obj.J, 'mean', Hm0_bins.centers, Te_bins.centers);

            plot_matrix(M,'Wave Energy Flux Matrix',"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_power_performance_workflow(testCase)
            filename = 'Capture Length Matrix mean.png';
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
            assertTrue(testCase,isfield(x,"mean"));
            assertTrue(testCase,isfield(x,"min"));
            assertTrue(testCase,isfield(x,"max"));
            assertTrue(testCase,isfield(x,"std"));
            assertTrue(testCase,isfield(x,"median"));
            assertTrue(testCase,isfield(x,"count"));
            assertTrue(testCase,isfield(x,"sum"));
            assertTrue(testCase,isfield(x,"freq"));
            assertEqual(testCase,y,401239.4822345051, 'RelTol',0.00001);

        end

    end

end
