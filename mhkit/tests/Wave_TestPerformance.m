classdef Wave_TestPerformance < matlab.unittest.TestCase

    methods (Test)

        function test_capture_length(testCase)
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);

            L = capture_length(Obj.P, Obj.J);
%             L_stats = statistics(L)
%             assertAlmostEqual(testCase,L_stats('mean'), 0.6676, 3);
        end

        function test_capture_length_matrix(testCase)
            seednum = 123;
            rng(seednum);
            Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);
            Obj.Hm0 = raylrnd(4, [1,100000]);
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            
            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te, L, 'std', Obj.Hm0_bins, Obj.Te_bins);
            assertEqual(testCase,size(LM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(LM.values))), 131);
        end

        function test_wave_energy_flux_matrix(testCase)
            seednum = 123;
            rng(seednum);
            Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);
            Obj.Hm0 = raylrnd(4, [1,100000]);
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            
            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            assertEqual(testCase,size(JM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(JM.values))), 131);
        end

        function test_power_matrix(testCase)
            seednum = 123;
            rng(seednum);
            Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);
            Obj.Hm0 = raylrnd(4, [1,100000]);
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            
            L = capture_length(Obj.P, Obj.J);
            LM = capture_length_matrix(Obj.Hm0, Obj.Te,L, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            JM = wave_energy_flux_matrix(Obj.Hm0, Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);
            PM = power_matrix(LM, JM);
            assertEqual(testCase,size(PM.values), [38 9]);
            assertEqual(testCase,sum(sum(isnan(PM.values))), 131);
        end

        function test_mean_annual_energy_production(testCase)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.TestCase
            import matlab.unittest.constraints.AbsoluteTolerance
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);
            
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
            Obj.Te = normrnd(4.5, 0.8, [1,100000]);
            Obj.P = normrnd(200, 40, [1,100000]);
            Obj.J = normrnd(300, 10, [1,100000]);
            Obj.Hm0 = raylrnd(4, [1,100000]);
            Obj.Hm0_bins = 0:0.5:18.5;
            Obj.Te_bins = 0:1:8;
            M = wave_energy_flux_matrix(Obj.Hm0,Obj.Te,Obj.J, 'mean', Obj.Hm0_bins, Obj.Te_bins);

            plot_matrix(M,'Wave Energy Flux Matrix',"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
    end
end
