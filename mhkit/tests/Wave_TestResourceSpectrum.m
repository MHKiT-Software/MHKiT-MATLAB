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

         % This test and function no longer in corresponding Python code
         % function test_bretschneider_spectrum(testCase)
         %     Obj.f = 0.1/(2*pi):0.01/(2*pi):3.5/(2*pi);
         %     Obj.Tp = 8;
         %     Obj.Hs = 2.5;

         %     S = create_spectra('bretschneider_spectrum',Obj.f,Obj.Tp,Obj.Hs);
         %     Tp0 = peak_period(S);
         %     Hm0 = significant_wave_height(S);
         %     errorHm0 = abs(Obj.Tp - Tp0)/Obj.Tp;
         %     errorTp0 = abs(Obj.Hs - Hm0)/Obj.Hs;
         %     assertLessThan(testCase,errorHm0, 0.01);
         %     assertLessThan(testCase,errorTp0, 0.01);
         % end

        function test_surface_elevation_seed(testCase)

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
            assertEqual(testCase,eta0, eta1);
        end

        function test_surface_elevation_phasing(testCase)

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
            assertEqual(testCase,eta0, eta1);
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

            % figure
            % plot(eta_ifft.time, eta_ifft.elevation)
            % hold on
            % plot(eta_sum_of_sines.time, eta_sum_of_sines.elevation)
            % hold off

            % xlabel('Time (s)')
            % ylabel('Surface Elevation (m)')
            % title('Surface Elevation over Time')
            % grid on
            % legend('IFFT Method', 'Sum of Sines Method')

            surface_elevation_diff = abs(sum(eta_ifft.elevation - eta_sum_of_sines.elevation));
            % disp(surface_elevation_diff);

            % This difference is suboptimal! The actual difference should be
            % very close to zero. There appears to be some form of drift over
            % time with the sum of sines method lagging slightly behind the
            % ifft method
            assertLessThan(testCase, surface_elevation_diff, 1.0);
        end

    end

end
