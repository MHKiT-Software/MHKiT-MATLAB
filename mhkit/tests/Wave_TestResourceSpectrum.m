classdef Wave_TestResourceSpectrum < matlab.unittest.TestCase

    methods (Test)

        function test_pierson_moskowitz_spectrum(testCase)
            Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;

            S = pierson_moskowitz_spectrum(Obj.f,Obj.Tp,Obj.Hs);
            Tp0 = peak_period(S);
            error = abs(Obj.Tp - Tp0)/Obj.Tp;
            assertLessThan(testCase,error, 0.01);
        end

        function test_surface_elevation_seed(testCase)
            assumeFail(testCase, "Per @simmsa: This test does not seem valid. A random seed and a defined seed should not have the same result")
            Obj.f = 0.0:0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.01/(2*pi);
            Trep = 1/df;
            Obj.t = 0:0.05:Trep;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            seednum = 123;
            eta0 = surface_elevation(S, Obj.t);
            eta1 = surface_elevation(S, Obj.t,"seed",seednum);
            assertEqual(testCase,eta0.elevation, eta1.elevation);
        end

        function test_surface_elevation_phasing(testCase)
            assumeFail(testCase, "Per @simmsa: This test does not account for randomness of the seed")
            Obj.f = 0.0:0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.01/(2*pi);
            Trep = 1/df;
            Obj.t = 0:0.05:Trep;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            eta0 = surface_elevation(S, Obj.t);
            seednum = 123;
            rng(seednum);
            phases = rand(size(S.spectrum))*2*pi;
            eta1 = surface_elevation(S, Obj.t,"phases",phases);
            assertEqual(testCase,eta0.elevation, eta1.elevation);
        end

        function test_surface_elevation_moments(testCase)

            Obj.f = 0.0:0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.01/(2*pi);
            Trep = 1/df;
            Obj.t = 0:0.05:Trep;
            dt = Obj.t(2)-Obj.t(1);

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            wave_elevation = surface_elevation(S, Obj.t);
            Sn = elevation_spectrum(wave_elevation.elevation, 1/dt,length(wave_elevation.elevation),Obj.t,"window","boxcar","detrend",false,"noverlap",0);
            m0 = frequency_moment(S,0);
            m0n = frequency_moment(Sn,0);
            errorm0 = abs((m0 - m0n)/m0);
            assertLessThan(testCase,errorm0, 0.01);
            m1 = frequency_moment(S,1);
            m1n = frequency_moment(Sn,1);
            errorm1 = abs((m1 - m1n)/m1);
            assertLessThan(testCase,errorm1, 0.01);
        end
%
%         function test_surface_elevation_rmse(testCase)
%             Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
%             Obj.Tp = 8;
%             Obj.Hs = 2.5;
%             df = 0.01/(2*pi);
%             Trep = 1/df;
%             Obj.t = 0:0.05:Trep;
%             import matlab.unittest.qualifications.Assertable
%
%             S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
%             wave_elevation = surface_elevation(S, Obj.t);
%             Sn = elevation_spectrum(wave_elevation.elevation,1/df,length(wave_elevation.elevation),Obj.t,"window","boxcar","detrend",false,"noverlap",0);
% %             fSn = interp1(Sn.frequency,Sn.spectrum,0:Sn.sample_rate:Sn.nnft);
% %             rmse = (S - fSn(Sn.frequency))^2;
% %             rmse_sum = (sum(rmse)/length(rmse))^0.5;
% %             assertLessThan(testCase,rmse_sum, 0.02);
%         end

        function test_jonswap_spectrum(testCase)
            Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            Hm0 = significant_wave_height(S);
            Tp0 = peak_period(S);
            errorHm0 = abs(Obj.Tp - Tp0)/Obj.Tp;
            errorTp0 = abs(Obj.Hs - Hm0)/Obj.Hs;
            assertLessThan(testCase,errorHm0, 0.01);
            assertLessThan(testCase,errorTp0, 0.01);
        end

        function test_jonswap_spectrum_gamma(testCase)
            Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            Obj.gamma = 2.0;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs, Obj.gamma);
            Hm0 = significant_wave_height(S);
            Tp0 = peak_period(S);
            errorHm0 = abs(Obj.Tp - Tp0)/Obj.Tp;
            errorTp0 = abs(Obj.Hs - Hm0)/Obj.Hs;
            assertLessThan(testCase,errorHm0, 0.01);
            assertLessThan(testCase,errorTp0, 0.01);
        end

        function test_plot_spectrum(testCase)
            Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;

            filename = 'wave_plot_matrix.png';
            if isfile(filename)
                delete(filename);
            end

            S = pierson_moskowitz_spectrum(Obj.f,Obj.Tp,Obj.Hs);

            plot_spectrum(S,"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function testSurfaceElevationMethod(testCase)
            Trep = 600;
            df = 1 / Trep;
            f = 0:df:1;
            Hs = 2.5;
            Tp = 8;
            t = 0:0.05:Trep;

            S = pierson_moskowitz_spectrum(f, Tp, Hs);

            eta_ifft = surface_elevation(S, t, "seed", 1, "method", "ifft");
            eta_sum_of_sines = surface_elevation(S, t, "seed", 1, "method", "sum_of_sines");

            surface_elevation_diff = mean(eta_ifft.elevation - eta_sum_of_sines.elevation);

            assertLessThan(testCase, surface_elevation_diff, 0.01);
        end

    end

end
