classdef Power_TestDevice < matlab.unittest.TestCase

    methods (Test)

        function test_harmonics_sine_wave(testCase)
            time = [0;1;2;3];
            t = 600;
            fs = 1000;
            samples = linspace(0,t,fs*t);
            frequency = 60;
            harmonics_int = 0:5:60*51;
            harmonics_int(end) = [];
            harmonics_vals = zeros(size(harmonics_int));
            harmonics_vals(13)= 1.0;  %setting 60th harmonic to amplitude of the signal
            current = sin(2*pi*frequency*samples);

            x = struct('current',current','time',samples');

            harmonic = harmonics(x, 1000, frequency);

            assertEqual(testCase,harmonic.amplitude', harmonics_vals,'AbsTol',0.1)
        end

        function test_harmonic_subgroup_sine_wave(testCase)
            frequency = 60;
            harmonics_int = 0:5:60*51;
            harmonics_int(end) = [];
            harmonics_vals = zeros(size(harmonics_int));
            harmonics_vals(13)= 1.0;  %setting 60th harmonic to amplitude of the signal
            harmonic_groups = harmonics_vals(1:12:end); %harmonic groups should be equal to every 12th harmonic in this idealized example

            harmonics = struct('amplitude',harmonics_vals','harmonic',harmonics_int);

            hsg = harmonic_subgroups(harmonics,frequency);

            assertEqual(testCase,hsg.amplitude', harmonic_groups,'AbsTol',0.1)
        end

        function test_TCHD_sine_wave(testCase)
            thcd = 0.0; %Since this is an idealized sin wave, there should be no distortion
            frequency = 60;
            harmonics_int = 0:5:60*60;
            harmonics_vals = zeros(size(harmonics_int));
            harmonics_vals(13)= 1.0;  %setting 60th harmonic to amplitude of the signal
            harmonics = struct('amplitude',harmonics_vals','harmonic',harmonics_int);

            hsg = harmonic_subgroups(harmonics,frequency);

            TCHD = total_harmonic_current_distortion(hsg,18.8); % had to just put a random rated current in here
            assertEqual(testCase,TCHD(1),thcd)
        end

        function test_interharmonics_sine_wave(testCase)
            frequency = 60;
            harmonics_int = 0:5:60*51;
            harmonics_int(end) = [];
            harmonics_vals = zeros(size(harmonics_int));
            harmonics_vals(13)= 1.0;  %setting 60th harmonic to amplitude of the signal
            harmonic_groups = harmonics_vals(1:12:end);
            interharmonic = zeros(size(harmonic_groups)); %since this is an idealized sin wave, the interharmonics should be zero
            harmonics = struct('amplitude',harmonics_vals','harmonic',harmonics_int);

            inter_harmonics = interharmonics(harmonics,frequency);

            % for i,j in zip(inter_harmonics.values, interharmonic):
            assertEqual(testCase,inter_harmonics.amplitude', interharmonic,'AbsTol',0.1)
        end

        function test_instfreq(testCase)
            frequency = 60;
            t = 600;
            fs = 1000;
            samples = 0:0.001:600;%linspace(0,t,fs*t)
            samples(end)=[];
            voltage = sin(2*pi*frequency*samples);
            si = size(voltage);
            inst_freq = ones(si(2)-1,si(1))*60;

            um = struct('voltage',voltage','time',samples);
            freq = instantaneous_frequency(um);

            assertEqual(testCase, freq.frequency, inst_freq,'AbsTol',0.1)
        end

        function test_dc_power(testCase)
            current_data = [1,2,3;4,5,6;7,8,9;10,11,12];
            voltage_data = [1,5,9;2,6,10;3,7,11;4,8,12];
            voltage.voltage = voltage_data;
            voltage.time = [0;1;2;3];
            current.current = current_data;
            current.time = [0;1;2;3];

            P = dc_power(voltage,current);
            assertEqual(testCase,sum(P.gross),sum(voltage.voltage.*current.current,'all'))
        end

        function test_ac_power_three_phase(testCase)
            current_data = [1,2,3;4,5,6;7,8,9;10,11,12];
            voltage_data = [1,5,9;2,6,10;3,7,11;4,8,12];
            voltage.voltage = voltage_data;
            voltage.time = [0;1;2;3];
            current.current = current_data;
            current.time = [0;1;2;3];

            P1 = ac_power_three_phase(voltage, current, 1, false);
            P1b = ac_power_three_phase(voltage, current, 0.5, false);
            P2 = ac_power_three_phase(voltage, current,1, true);
            P2b = ac_power_three_phase(voltage, current, 0.5, true);

            assertEqual(testCase,sum(P1.power), 584)
            assertEqual(testCase,sum(P1b.power), 584/2)
            assertEqual(testCase,sum(P2.power), 1011.518,'AbsTol', 0.01)
            assertEqual(testCase,sum(P2b.power), 1011.518/2,'AbsTol', 0.01)
        end

        function test_gen_test_data(testCase)
            % Sr, Un, In, SCR, fg, & fs
            Sr = 3e6; Un=12e3; In=144; SCR=20; fg=60; fs=50e3;fv=0.5;
            % 1. opt=0, pure sine wave with alpha0=pi/6
            opt = 0;fm = 20;
            DeltaI_I = [0 0 0 0];% pure sine wave
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,10);
            u_m0 = readmatrix("../../examples/data/power/testdata/sinewave-pi6_um.txt");
            i_m0 = readmatrix("../../examples/data/power/testdata/sinewave-pi6_im.txt");
            testCase.verifyTrue(max(abs( ...
                (i_m.data(1:100,1)-i_m0(:,1))./i_m0(:,1) ...
                ))<1e-10,opt);
            testCase.verifyTrue(max(abs( ...
                (u_m.data(1:100)-u_m0)./u_m0)) ...
                <1e-10,opt);

            % 2. opt=randi([2,5]) tests generated according to IECTS
            opt = randi([2,5],1);
            if opt==2
                SCR = 50; fm = 33.3;
                % TableB.3,fg=60, SCR=50
                DeltaI_I = [8.040 9.376 11.479 12.844]; % fm=33.3
            else
                SCR = 20; fm = 25;
                % TableB.2, fg=60, SCR=20
                DeltaI_I = [4.763 5.726 7.640 9.488];% fm=25
                %DeltaI_I = [3.212 3.958 5.644 7.711];% fm=20
            end
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,10);
            i_m0 = readmatrix( ...
                sprintf('../../examples/data/power/testdata/B.3.%i_im.txt',opt));
            u_m0 = readmatrix( ...
                sprintf('../../examples/data/power/testdata/B.3.%i_um.txt',opt));
            testCase.verifyTrue(max(abs( ...
                (i_m.data(1:100,1)-i_m0(:,1))./i_m0(:,1) ...
                ))<1e-10,opt);
            testCase.verifyTrue(max(abs((u_m.data(1:100)-u_m0)./u_m0))<1e-10,...
                opt);
        end

        function test_flicker_ufic_workflow(testCase)
            % u_m, i_m, Sr, Un, In, SCR, fg, & fs
            Sr = 3e6; Un=12e3; In=144; fg=60; fs=50e3;fv=0.5;
            %1. opt = 0, pure sine waves:
            opt = 0; idx = randi([1,4],1);
            fm=20; SCR=20; DeltaI_I = [0 0 0 0];% pure sine wave
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,10);
            method = 'ZCD'; methodopts = {};
            out = flicker_ufic_workflow( ...
                Sr,Un,SCR,fg,u_m,i_m,method,methodopts);
            freq0=readmatrix("../../examples/data/power/testdata/sinewave-pi6_SCR20_freq.txt");
            alpha_m0 =readmatrix("../../examples/data/power/testdata/sinewave-pi6_SCR20_alpham.txt");
            u00=readmatrix("../../examples/data/power/testdata/sinewave-pi6_SCR20_u0.txt");
            u_fic0=readmatrix("../../examples/data/power/testdata/sinewave-pi6_SCR20_ufic.txt");
            testCase.verifyTrue(max(abs( ...
                (out.freq.data(1:100)-freq0)./freq0))<1e-10,string(opt));
            testCase.verifyTrue(max(abs( ...
                (out.alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(max(abs( ...
                (out.u0(1:100)-u00)./u00))<1e-10,string(opt));
            testCase.verifyTrue(max(abs(( ...
                out.u_fic(1:100,idx)-u_fic0(:,idx))./u_fic0(:,idx)))<1e-10, ...
                string(opt));
            % 2. opt=randi([2,5]) tests generated according to IECTS
            opt = randi([2,5],1); idx = randi([1,4],1);
            if opt==2
                SCR = 50; fm = 33.3;
                % TableB.3,fg=60,SCR=50
                DeltaI_I = [8.040 9.376 11.479 12.844]; % fm=33.3
            else
                SCR = 20; fm = 25;
                % TableB.2,fg=60,SCR=20
                DeltaI_I = [4.763 5.726 7.640 9.488];% fm=25
                %DeltaI_I = [3.212 3.958 5.644 7.711];% fm=20
            end
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,10);
            if opt==3
                % B.3.3 Distorted um with multiple zero crossings
                method = 'stft';methodopts = {...
                    'Window',rectwin(int32(50000)),'OverlapLength',25000,...
                    'FFTLength',50e3,'FrequencyRange','onesided'};
            else
                method = 'ZCD'; methodopts = {};
            end
            out = flicker_ufic_workflow(...
                Sr,Un,SCR,fg,u_m,i_m,method,methodopts);
            freq0    = readmatrix( ...
                sprintf('../../examples/data/power/testdata/B.3.%i_freq.txt',opt));
            alpha_m0 = readmatrix( ...
                sprintf('../../examples/data/power/testdata/B.3.%i_alpham.txt',opt));
            u00      = readmatrix( ...
                sprintf('../../examples/data/power/testdata/B.3.%i_u0.txt',opt));
            u_fic0   = readmatrix(...
                sprintf('../../examples/data/power/testdata/B.3.%i_ufic.txt',opt));
            testCase.verifyTrue(max(abs( ...
                (out.freq.data(1:100)-freq0)./freq0))<1e-10, string(opt));
            testCase.verifyTrue(max(abs( ...
                (out.alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10, string(opt));
            testCase.verifyTrue(max(abs( ...
                (out.u0(1:100)-u00)./u00))<1e-10, string(opt));
            testCase.verifyTrue(max(abs(( ...
                out.u_fic(1:100,idx)-u_fic0(:,idx))./u_fic0(:,idx)))<1e-10, ...
                string(opt));
        end

        function test_calc_Rfic_Lfic(testCase)
            Sr = 3e6; Un=12e3; fg=60;
            % 1. SCR = 20:
            SCR = 20;
            [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
            Rfic0=readmatrix("../../examples/data/power/testdata/FicGrid_SCR20_fg60_Rfic.txt");
            Lfic0=readmatrix("../../examples/data/power/testdata/FicGrid_SCR20_fg60_Lfic.txt");
            testCase.verifyTrue(max(abs((Rfic-Rfic0)./Rfic0))<1e-10,string(SCR));
            testCase.verifyTrue(max(abs((Lfic-Lfic0)./Lfic0))<1e-10,string(SCR));
            % 2. SCR = 50:
            SCR = 50;
            [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
            Rfic0=readmatrix("../../examples/data/power/testdata/FicGrid_SCR50_fg60_Rfic.txt");
            Lfic0=readmatrix("../../examples/data/power/testdata/FicGrid_SCR50_fg60_Lfic.txt");
            testCase.verifyTrue(max(abs((Rfic-Rfic0)./Rfic0))<1e-10,string(SCR));
            testCase.verifyTrue(max(abs((Lfic-Lfic0)./Lfic0))<1e-10,string(SCR));
        end

        function test_calc_flicker_coefficient(testCase)
            P_stfic = 0.1016;
            S_kfic = 6e7;
            Sr = 3e6;
            coef_flicker = calc_flicker_coefficient(P_stfic,S_kfic,Sr);
            assertEqual(testCase,coef_flicker, 2.032);
        end

        function test_calc_shortterm_flicker_severity(testCase)
            % prep input:
            P = struct();P.p0p1=0.4073;
            P.p0p7=0.3242; P.p1=0.3236;P.p1p5 = 0.3230;
            P.p2p2=0.2620; P.p3=0.2615;P.p4=0.2610;
            P.p6=0.1628; P.p8=0.1625;P.p10=0.1619;P.p13=0.1613;P.p17=0.161;
            P.p30=0.0308;P.p50=0.0301;P.p80=0.0294;
            % calculate Pst
            Pst = calc_shortterm_flicker_severity(P);
            testCase.verifyTrue(abs((Pst-0.307731652580621)./Pst)<1e-6,'Pst');
        end

    end

end

