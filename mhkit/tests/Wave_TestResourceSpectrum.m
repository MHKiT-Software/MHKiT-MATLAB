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

        function test_pierson_moskowitz_spectrum_zero_freq(testCase)
            % f=0 should always evaluate to exactly 0 spectral density.
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L146-L155
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.1;
            f_zero = 0:df:1-df;
            f_nonzero = df:df:1-df;

            S_zero = pierson_moskowitz_spectrum(f_zero, Obj.Tp, Obj.Hs);
            S_nonzero = pierson_moskowitz_spectrum(f_nonzero, Obj.Tp, Obj.Hs);

            assertEqual(testCase, S_zero.spectrum(1), 0.0);
            assertGreaterThan(testCase, S_nonzero.spectrum(1), 0.0);
        end

        function test_jonswap_spectrum_zero_freq(testCase)
            % f=0 should always evaluate to exactly 0 spectral density.
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L204-L213
            Obj.Tp = 8;
            Obj.Hs = 2.5;
            df = 0.1;
            f_zero = 0:df:1-df;
            f_nonzero = df:df:1-df;

            S_zero = jonswap_spectrum(f_zero, Obj.Tp, Obj.Hs);
            S_nonzero = jonswap_spectrum(f_nonzero, Obj.Tp, Obj.Hs);

            assertEqual(testCase, S_zero.spectrum(1), 0.0);
            assertGreaterThan(testCase, S_nonzero.spectrum(1), 0.0);
        end

        function test_surface_elevation_uses_sum_of_sines_when_input_frequency_index_does_not_have_zero(testCase)
            % ifft requires f(1)==0; without it, the default method must
            % silently fall back to sum_of_sines and match it exactly.
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/tests/wave/test_resource_spectrum.py#L216-L227
            f = linspace(1/30, 1/2, 32);
            Hs = 2.5;
            Tp = 8;
            Trep = 600;
            t = 0:0.05:Trep-0.05; % matches np.arange(0, Trep, 0.05), which excludes the endpoint

            S = jonswap_spectrum(f, Tp, Hs);

            eta_default = surface_elevation(S, t, "seed", 1);
            eta_sos = surface_elevation(S, t, "seed", 1, "method", "sum_of_sines");

            assertTrue(testCase, isfield(eta_default, 'elevation'));
            assertEqual(testCase, eta_default.elevation, eta_sos.elevation, 'AbsTol', 1e-6);
        end

        function test_surface_elevation_warn_user_if_zero_frequency_not_defined_and_using_ifft(testCase)
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/tests/wave/test_resource_spectrum.py#L229-L236
            f = linspace(1/30, 1/2, 32);
            Hs = 2.5;
            Tp = 8;
            Trep = 600;
            t = 0:0.05:Trep-0.05; % matches np.arange(0, Trep, 0.05), which excludes the endpoint

            S = jonswap_spectrum(f, Tp, Hs);

            testCase.verifyWarning(@() surface_elevation(S, t, "seed", 1, "method", "ifft"), ...
                'MHKiT:surface_elevation:MethodFallback');
        end

        function test_surface_elevation_uses_ifft_when_input_frequency_index_has_zero(testCase)
            % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/tests/wave/test_resource_spectrum.py#L238-L244
            Trep = 600;
            df = 1 / Trep;
            f = 0:df:1-df; % matches np.arange(0, 1, df), which excludes the endpoint
            Hs = 2.5;
            Tp = 8;
            t = 0:0.05:Trep-0.05; % matches np.arange(0, Trep, 0.05), which excludes the endpoint

            S = jonswap_spectrum(f, Tp, Hs);

            eta_default = surface_elevation(S, t, "seed", 1);
            eta_ifft = surface_elevation(S, t, "seed", 1, "method", "ifft");

            assertEqual(testCase, eta_default.elevation, eta_ifft.elevation);
        end
    end

end
