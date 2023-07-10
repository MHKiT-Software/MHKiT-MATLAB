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
%             end
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

%             for i,j in zip(inter_harmonics.values, interharmonic):
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
            % u_m, i_m, Sr, Un, In, SCR, fg, & fs
            Sr = 3e6; Un=12e3; In=144; SCR=20; fg=60; fs=50e3;
            % 1. opt=0, pure sine wave
            opt = 0;fm = 20;
            DeltaI_I = [0 0 0 0];% pure sine wave
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,10);
            u_m0 = readmatrix("testdata/sinewave-pi6_um.txt");
            i_m0 = readmatrix("testdata/sinewave-pi6_im.txt");
            testCase.verifyTrue(max(abs( ...
                (i_m.data(1:100,1)-i_m0(:,1))./i_m0(:,1) ...
                ))<1e-10,opt);
            testCase.verifyTrue(max(abs((u_m.data(1:100)-u_m0)./u_m0))<1e-10,...
                opt);
            
            % 2. opt=4, Distorted um with inter-harmonic modulation
            opt = 4; SCR=20; fm = 25;
            DeltaI_I = [4.763 5.726 7.640 9.488];% table B.2
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,10);
            i_m0 = readmatrix("testdata/B.3.4_fm25_SCR20_im.txt");
            u_m0 = readmatrix("testdata/B.3.4_fm25_SCR20_fv0.5_um.txt");
            testCase.verifyTrue(max(abs( ...
                (i_m.data(1:100,1)-i_m0(:,1))./i_m0(:,1) ...
                ))<1e-10,opt);
            testCase.verifyTrue(max(abs((u_m.data(1:100)-u_m0)./u_m0))<1e-10,...
                opt);
            % 2. opt=5, B.3.5 slow frequency changes
            opt = 5; SCR=20; fm = 25;
            %DeltaI_I = [3.212 3.958 5.644 7.711]; % fm=20
            DeltaI_I = [4.763 5.726 7.640 9.488];% table B.2 fm=25
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,10);
            i_m0 = readmatrix("testdata/B.3.5_fm25_SCR20_im.txt");
            u_m0 = readmatrix("testdata/B.3.5_fm25_SCR20_um.txt");
            testCase.verifyTrue(max(abs( ...
                (i_m.data(1:100,2)-i_m0(:,2))./i_m0(:,2) ...
                ))<1e-10,opt);
            testCase.verifyTrue(max(abs((u_m.data(1:100)-u_m0)./u_m0))<1e-10,...
                opt);
        end

        function test_flicker_assessment(testCase)
            %1 opt = 0, pure sine waves:
            % u_m, i_m, Sr, Un, In, SCR, fg, & fs
            Sr = 3e6; Un=12e3; In=144; fg=60; fs=50e3;
            fm = 20;SCR=20;opt = 0;
            DeltaI_I = [0 0 0 0];% pure sine wave
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,650);
            method = 'ZCD'; opts = {};
            [~,freq,alpha_m,u0,u_fic] = ...
                flicker_assessment(Sr,Un,SCR,fg,u_m,i_m,method,opts);
            freq0=readmatrix("testdata/sinewave-pi6_SCR20_freq.txt");
            alpha_m0 =readmatrix("testdata/sinewave-pi6_SCR20_alpham.txt");
            u00=readmatrix("testdata/sinewave-pi6_SCR20_u0.txt");
            u_fic0=readmatrix("testdata/sinewave-pi6_SCR20_ufic.txt");
            testCase.verifyTrue(...
                max(abs(freq.data(1:100)-freq0)./freq0(:,1))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u0(1:100)-u00)./u00(:,1)))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u_fic(1:100,2)-u_fic0(:,2))./u_fic0(:,2)))<1e-10,string(opt));
            % 2 opt = 2, B.3.2
            fm = 33.3;SCR=50;opt = 2;
            % TableB.3,fg=60, SCR=50
            DeltaI_I = [8.040 9.376 11.479 12.844]; % fm=33.3
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,650);
            method = 'ZCD'; opts = {};
            [~,freq,alpha_m,u0,u_fic] = ...
                flicker_assessment(Sr,Un,SCR,fg,u_m,i_m,method,opts);
            freq0=readmatrix(    "testdata/B.3.2_fm33.3_SCR50_freq.txt");
            alpha_m0 =readmatrix("testdata/B.3.2_fm33.3_SCR50_alpham.txt");
            u00=readmatrix(      "testdata/B.3.2_fm33.3_SCR50_u0.txt");
            u_fic0=readmatrix(   "testdata/B.3.2_fm33.3_SCR50_ufic.txt");
            testCase.verifyTrue(...
                max(abs(freq.data(1:100)-freq0)./freq0(:,1))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u0(1:100)-u00)./u00(:,1)))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u_fic(1:100,1)-u_fic0(:,1))./u_fic0(:,1)))<1e-10, string(opt));
            % 3 opt = 3, B.3.3 Distorted um with multiple zero crossings
            fm = 25;SCR=20;opt = 3;
            % TableB.2,fg=60, SCR=20
            DeltaI_I = [4.763 5.726 7.640 9.488]; % fm=25
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,650);
            method = 'STFT'; 
            opts = {'Window',rectwin(50000),'OverlapLength',25000,...
                'FFTLength',50e3,'FrequencyRange','onesided'};
            [~,freq,alpha_m,u0,u_fic] = ...
                flicker_assessment(Sr,Un,SCR,fg,u_m,i_m,method,opts);
            freq0=readmatrix(    "testdata/B.3.3_fm25_SCR20_freq.txt");
            alpha_m0 =readmatrix("testdata/B.3.3_fm25_SCR20_alpham.txt");
            u00=readmatrix(      "testdata/B.3.3_fm25_SCR20_u0.txt");
            u_fic0=readmatrix(   "testdata/B.3.3_fm25_SCR20_ufic.txt");
            testCase.verifyTrue(...
                max(abs(freq.data(1:100)-freq0)./freq0(:,1))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u0(1:100)-u00)./u00(:,1)))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u_fic(1:100,3)-u_fic0(:,3))./u_fic0(:,3)))<1e-10, string(opt));
            % 4 opt = 4, B.3.4 Distorted um with inter-harmonic modulation
            fm = 25;SCR=20;opt = 4;
            % TableB.2,fg=60, SCR=20
            DeltaI_I = [4.763 5.726 7.640 9.488]; % fm=25
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,650);
            method = 'zcd'; 
            [~,freq,alpha_m,u0,u_fic] = ...
                flicker_assessment(Sr,Un,SCR,fg,u_m,i_m,method,opts);
            freq0=readmatrix(    "testdata/B.3.4_fm25_SCR20_freq.txt");
            alpha_m0 =readmatrix("testdata/B.3.4_fm25_SCR20_alpham.txt");
            u00=readmatrix(      "testdata/B.3.4_fm25_SCR20_u0.txt");
            u_fic0=readmatrix(   "testdata/B.3.4_fm25_SCR20_ufic.txt");
            testCase.verifyTrue(...
                max(abs(freq.data(1:100)-freq0)./freq0(:,1))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u0(1:100)-u00)./u00(:,1)))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u_fic(1:100,3)-u_fic0(:,3))./u_fic0(:,3)))<1e-10, string(opt));
            % 5 opt = 5, B.3.5 slow frequency changes
            fm = 25;SCR=20;opt = 5;
            % TableB.2,fg=60, SCR=20
            DeltaI_I = [4.763 5.726 7.640 9.488]; % fm=25
            [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,650);
            method = 'zcd'; 
            [~,freq,alpha_m,u0,u_fic] = ...
                flicker_assessment(Sr,Un,SCR,fg,u_m,i_m,method,opts);
            freq0=readmatrix(    "testdata/B.3.5_fm25_SCR20_freq.txt");
            alpha_m0 =readmatrix("testdata/B.3.5_fm25_SCR20_alpham.txt");
            u00=readmatrix(      "testdata/B.3.5_fm25_SCR20_u0.txt");
            u_fic0=readmatrix(   "testdata/B.3.5_fm25_SCR20_ufic.txt");
            testCase.verifyTrue(...
                max(abs(freq.data(1:100)-freq0)./freq0(:,1))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((alpha_m(1:100)-alpha_m0)./alpha_m0))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u0(1:100)-u00)./u00(:,1)))<1e-10,string(opt));
            testCase.verifyTrue(...
                max(abs((u_fic(1:100,4)-u_fic0(:,4))./u_fic0(:,4)))<1e-10, string(opt));
        end
        function test_calc_Rfic_Lfic(testCase)
            Sr = 3e6; Un=12e3; fg=60; 
            % 1. SCR = 20:
            SCR = 20;
            [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
            Rfic0=readmatrix("testdata/FicGrid_SCR20_fg60_Rfic.txt");
            Lfic0=readmatrix("testdata/FicGrid_SCR20_fg60_Lfic.txt");
            testCase.verifyTrue(max(abs((Rfic-Rfic0)./Rfic0))<1e-10,string(SCR));
            testCase.verifyTrue(max(abs((Lfic-Lfic0)./Lfic0))<1e-10,string(SCR));
            % 2. SCR = 50:
            SCR = 50;
            [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
            Rfic0=readmatrix("testdata/FicGrid_SCR50_fg60_Rfic.txt");
            Lfic0=readmatrix("testdata/FicGrid_SCR50_fg60_Lfic.txt");
            testCase.verifyTrue(max(abs((Rfic-Rfic0)./Rfic0))<1e-10,string(SCR));
            testCase.verifyTrue(max(abs((Lfic-Lfic0)./Lfic0))<1e-10,string(SCR));
        end
    end
end  




        