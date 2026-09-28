classdef Acoustics_TestMetrics < matlab.unittest.TestCase
    properties
        spsd
        spsd_60s
    end

    methods (TestClassSetup)
        function setupOnce(testCase)
            testdir = fileparts(mfilename('fullpath'));
            datadir = fullfile(testdir, '..', '..', 'examples', 'data', 'acoustics');
            file_name = fullfile(datadir, '6247.230204150508.wav');
            P = read_soundtrap(file_name, -177);
            testCase.spsd = sound_pressure_spectral_density(P, P.fs, 1);
            testCase.spsd_60s = sound_pressure_spectral_density(P, P.fs, 60);
        end
    end

    methods (Test)
        function test_spl(testCase)
            td_spl = sound_pressure_level(testCase.spsd, 10, 100000);
            td_spl10 = decidecade_sound_pressure_level(testCase.spsd, 10, 100000);
            td_spl3 = third_octave_sound_pressure_level(testCase.spsd, 10, 100000);

            cc = datetime([ ...
                "2023-02-04 15:05:08.000000"
                "2023-02-04 15:05:08.501683"
                "2023-02-04 15:05:09.003367"
                "2023-02-04 15:05:09.505051"
                "2023-02-04 15:05:10.006735"], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
            cd_spl_head = [98.1252, 98.3950, 98.6400, 98.3964, 97.6300];
            cd_spl_tail = [97.43068, 97.67616, 97.99467, 98.13005, 97.96256];

            cd_spl10_freq_head = [10.0, 12.589254, 15.848932, 19.952623, 25.118864];
            cd_spl10_head = [
                67.698427107731632	68.718504549486369	55.150239136134928	63.092171166170736	68.285316609698185;
                70.638117562215356	68.352358818755633	70.679291272194590	66.826965002510178	72.307235015440625;
                75.061041066120794	72.000060887346237	61.867764868621322	69.029564442104871	65.566121037342981;
                73.571073481923875	77.646357542405568	71.315874138184853	68.348127773527210	74.641385774888377;
                83.318133484432281	82.382487625370942	82.839583740048681	80.066178195328149	82.922058814906620
            ];
            cd_spl10_freq_tail = [19952.62315, 25118.864315, 31622.776602, 39810.717055, 50118.723363];
            cd_spl10_tail = [
                81.263857874463184	81.623440496905090	81.705362821962368	82.055563841890347	80.905483744451274;
                81.411568612914579	81.499648040350351	81.424401084480735	81.612600277281615	81.397403019626026;
                83.528128225083123	83.623441861246164	83.454856374486340	83.690034791494767	83.367940400050202;
                81.817594809639502	81.759206900301180	81.471336683096752	81.647659460075118	81.573892228707166;
                74.153354148096412	74.155760294865715	73.855762818092344	74.059185929338497	74.344841337560410
            ];
            cd_spl3_freq_head = [10.0, 12.59921, 15.874011, 20.0, 25.198421];
            cd_spl3_head = [
                67.698427107731632	68.718504549486369	55.150239136134928	63.092171166170736	68.285316609698185;
                70.638117562215356	68.352358818755633	70.679291272194590	66.826965002510178	72.307235015440625;
                75.061041066120794	72.000060887346237	61.867764868621322	69.029564442104871	65.566121037342981;
                73.571073481923875	77.646357542405568	71.315874138184853	68.348127773527210	74.641385774888377;
                83.318133484432281	82.382487625370942	82.839583740048681	80.066178195328149	82.922058814906620
            ];
            cd_spl3_freq_tail = [20480.0, 25803.183102, 32509.973544, 40960.0, 51606.366204];
            cd_spl3_tail = [
                81.134544889533359	81.479943775144307	81.531385268188657	81.881263944278700	80.703903519567945;
                81.677810946590085	81.775935532186026	81.694836721439614	81.865021419436431	81.690848660849610;
                83.968848694011569	83.970789087350326	83.796709880206919	84.065714823337430	83.760770230467969;
                80.663402401843385	80.686137257221461	80.384396424319945	80.522388118409964	80.532964378490249;
                72.079946592718045	72.278741310892528	71.948646690849614	71.850495761875621	72.240261721430514
            ];

            testCase.verifyEqual(td_spl.data(1:5), cd_spl_head, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl.data(end-4:end), cd_spl_tail, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl10.freq(1:5), cd_spl10_freq_head, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl10.data(1:5,1:5), cd_spl10_head, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl10.freq(end-4:end), cd_spl10_freq_tail, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl10.data(end-4:end,end-4:end), cd_spl10_tail, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl3.freq(1:5), cd_spl3_freq_head, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl3.data(1:5,1:5), cd_spl3_head, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl3.freq(end-4:end), cd_spl3_freq_tail, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_spl3.data(end-4:end,end-4:end), cd_spl3_tail, 'AbsTol', 1e-3);
            testCase.verifyEqual(posixtime(td_spl.time(1:5)), posixtime(cc'), 'AbsTol', 1e-5);
        end

        function test_nmfs_weighting(testCase)
            freq = testCase.spsd.freq;
            slc = 20:25; % MATLAB is 1-based

            [W_LF, E_LF] = nmfs_auditory_weighting(freq, 'LF');
            [W_HF, E_HF] = nmfs_auditory_weighting(freq, 'HF');
            [W_VHF, E_VHF] = nmfs_auditory_weighting(freq, 'VHF');
            [W_PW, E_PW] = nmfs_auditory_weighting(freq, 'PW');
            [W_OW, E_OW] = nmfs_auditory_weighting(freq, 'OW');

            cd_W_LF = [-18.241247, -17.827854, -17.434275, -17.058767, -16.699821, -16.3561];
            cd_E_LF = [195.36125, 194.94786, 194.55428, 194.17877, 193.81982, 193.4761];
            cd_W_HF = [-59.7284, -59.071625, -58.44541, -57.847057, -57.274178, -56.724693];
            cd_E_HF = [241.0484, 240.39163, 239.76541, 239.16705, 238.59418, 238.0447];
            cd_W_VHF = [-109.34241, -108.397385, -107.49632, -106.635315, -105.81097, -105.02029];
            cd_E_VHF = [270.2524, 269.30737, 268.4063, 267.54532, 266.72098, 265.9303];
            cd_W_PW = [-52.117348, -51.427025, -50.768852, -50.13999, -49.537937, -48.96051];
            cd_E_PW = [227.40735, 226.71703, 226.05885, 225.43, 224.82794, 224.25052];
            cd_W_OW = [-65.056496, -64.386955, -63.748577, -63.138584, -62.55456, -61.99438];
            cd_E_OW = [244.4265, 243.75696, 243.11858, 242.50858, 241.92456, 241.36438];

            testCase.verifyEqual(W_LF(slc), cd_W_LF', 'AbsTol', 1e-4);
            testCase.verifyEqual(W_HF(slc), cd_W_HF', 'AbsTol', 1e-4);
            testCase.verifyEqual(W_VHF(slc), cd_W_VHF', 'AbsTol', 1e-4);
            testCase.verifyEqual(W_PW(slc), cd_W_PW', 'AbsTol', 1e-4);
            testCase.verifyEqual(W_OW(slc), cd_W_OW', 'AbsTol', 1e-4);

            testCase.verifyEqual(E_LF(slc), cd_E_LF', 'AbsTol', 1e-4);
            testCase.verifyEqual(E_HF(slc), cd_E_HF', 'AbsTol', 1e-4);
            testCase.verifyEqual(E_VHF(slc), cd_E_VHF', 'AbsTol', 1e-4);
            testCase.verifyEqual(E_PW(slc), cd_E_PW', 'AbsTol', 1e-4);
            testCase.verifyEqual(E_OW(slc), cd_E_OW', 'AbsTol', 1e-4);
        end

        function test_sel(testCase)
            td_sel = sound_exposure_level(testCase.spsd_60s,[], 10, 100000);
            td_sel_lf = sound_exposure_level(testCase.spsd_60s,'LF', 10, 100000);
            td_sel_hf = sound_exposure_level(testCase.spsd_60s,'HF', 10, 100000);
            td_sel_vhf = sound_exposure_level(testCase.spsd_60s,'VHF', 10, 100000 );
            td_sel_pw = sound_exposure_level(testCase.spsd_60s,'PW', 10, 100000);
            td_sel_ow = sound_exposure_level(testCase.spsd_60s,'OW', 10, 100000);

            cc = datetime([ ...
                "2023-02-04 15:05:08.000000"
                "2023-02-04 15:05:45.500874"
                "2023-02-04 15:06:23.001750"
                "2023-02-04 15:07:00.502625"
                "2023-02-04 15:07:38.003500"], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
            cd_sel = [1.161832567214382e+02	1.177608180325696e+02	1.216988394734190e+02	1.322482028621648e+02	1.432811329170941e+02];
            cd_sel_lf = [1.123633263210942e+02	1.150560720648516e+02	1.201771340004273e+02	1.315286929637698e+02	1.427492968479198e+02];
            cd_sel_hf = [1.122226246068385e+02	1.143795287244090e+02	1.188811398942907e+02	1.293576472503867e+02	1.399411804663656e+02];
            cd_sel_vhf = [1.102333042133257e+02	1.111969207330859e+02	1.140075568256765e+02	1.228731328944110e+02	1.332000231806384e+02];
            cd_sel_pw = [1.122223727653538e+02	1.148390026535553e+02	1.198729727976822e+02	1.308165812210123e+02	1.416746474465080e+02];
            cd_sel_ow = [1.109456496065601e+02	1.133797705458589e+02	1.180640501702054e+02	1.285908810877421e+02	1.390843177149019e+02];

            testCase.verifyEqual(td_sel.data(1:5), cd_sel, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_sel_lf.data(1:5), cd_sel_lf, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_sel_hf.data(1:5), cd_sel_hf, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_sel_vhf.data(1:5), cd_sel_vhf, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_sel_pw.data(1:5), cd_sel_pw, 'AbsTol', 1e-3);
            testCase.verifyEqual(td_sel_ow.data(1:5), cd_sel_ow, 'AbsTol', 1e-3);
            testCase.verifyEqual(posixtime(td_sel.time(1:5)), posixtime(cc'), 'AbsTol', 1e-6);
        end

        function test_spl_vs_sel(testCase)
            td_spl = sound_pressure_level(testCase.spsd, 10, 100000);
            td_sel = sound_exposure_level(testCase.spsd, [], 10, 100000);
            testCase.verifyEqual(td_spl.data, td_sel.data, 'AbsTol', 1e-6);
        end
    end
end