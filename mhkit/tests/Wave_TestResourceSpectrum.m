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
            % An explicit seed must be reproducible across calls. MHKiT-Python
            % defaults to seed=None, so an unseeded call is not reproducible:
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L250
            Obj.f = 0.0:0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.01/(2*pi);
            Trep = 1/df;
            Obj.t = 0:0.05:Trep;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            seednum = 123;
            eta0 = surface_elevation(S, Obj.t, "seed", seednum);
            eta1 = surface_elevation(S, Obj.t, "seed", seednum);
            assertEqual(testCase,eta0.elevation, eta1.elevation);
        end

        function test_surface_elevation_phasing(testCase)
            % Explicit phases must match phases generated internally from the
            % same seed, following MHKiT-Python convention:
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L372-L380
            Obj.f = 0.0:0.01/(2*pi):3.5/(2*pi);
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.01/(2*pi);
            Trep = 1/df;
            Obj.t = 0:0.05:Trep;

            S = jonswap_spectrum(Obj.f, Obj.Tp, Obj.Hs);
            seednum = 123;
            eta0 = surface_elevation(S, Obj.t, "seed", seednum);
            rng(seednum);
            phases = rand(size(S.spectrum))*2*pi;
            eta1 = surface_elevation(S, Obj.t,"phases",phases);
            assertEqual(testCase,eta0.elevation, eta1.elevation);
        end

%

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
