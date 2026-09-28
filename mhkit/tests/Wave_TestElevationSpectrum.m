classdef Wave_TestElevationSpectrum < matlab.unittest.TestCase
    % Tests that depend on elevation_spectrum (pwelch), still Python-backed; native rewrite is part 2

    methods (Test)

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

    end

end
