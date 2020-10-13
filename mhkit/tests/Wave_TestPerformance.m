classdef Wave_TestPerformance < matlab.unittest.TestCase

    methods (Test)

        function test_capture_length(testCase)
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

        function test_capture_length_matrix(testCase)
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
            
            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te, L, 'std', Obj.Hm0_bins, Obj.Te_bins);
            
            assertEqual(testCase,size(LM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(LM.values))), 34, 'RelTol',0.1);
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
            assertEqual(testCase,sum(sum(isnan(JM.values))), 34, 'RelTol',0.1);
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
            
            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te,L, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            PM = power_matrix(LM, JM);
            assertEqual(testCase,size(PM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(PM.values))), 34, 'RelTol',0.1);
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
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            M = wave_energy_flux_matrix(Obj.Hm0,Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);

            plot_matrix(M,'Wave Energy Flux Matrix',"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
    end
end
