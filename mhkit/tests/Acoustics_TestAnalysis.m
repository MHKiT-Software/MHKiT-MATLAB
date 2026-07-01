classdef Acoustics_TestAnalysis < matlab.unittest.TestCase
    properties
        spsd
    end

    methods (TestClassSetup)
        function setupOnce(testCase)
            testdir = fileparts(mfilename('fullpath'));
            datadir = fullfile(testdir, '..', '..', 'examples', 'data', 'acoustics');
            file_name = fullfile(datadir, '6247.230204150508.wav');
            P = read_soundtrap(file_name, -177);
            testCase.spsd = sound_pressure_spectral_density(P, P.fs, 1);
        end
    end

    methods (Test)
        function test_sound_pressure_spectral_density(testCase)
            time = 0:0.1:9.9;
            data = sin(time);
            pressure.data = data;
            pressure.time = time;
            pressure.fs = 100;
            pressure.units = 'Pa';

            fs = 100;
            bin_length = 0.1;
            win_samples = bin_length * fs;

            spsd = sound_pressure_spectral_density(pressure, fs, bin_length);

            testCase.assertTrue(isstruct(spsd));
            testCase.assertTrue(isfield(spsd, 'freq'));
            testCase.assertTrue(isfield(spsd, 'time'));
            testCase.assertEqual(spsd.units, "Pa^2/Hz");
            testCase.assertEqual(spsd.bin_length, bin_length);
            testCase.assertEqual(spsd.n_fft, bin_length*fs);
            testCase.assertEqual(spsd.name, "Mean Square Sound Pressure Spectral Density");
            
            overlap = 0.5; % 50% overlap
            expected_segments = floor(length(pressure.data) / (win_samples*(1-overlap)));
            testCase.assertEqual(size(spsd.data, 2), expected_segments - 1);

            % Test with rms = false
            spsd_no_rms = sound_pressure_spectral_density(pressure, fs, bin_length, 0, false);
            testCase.assertEqual(spsd_no_rms.name, "Sound Pressure Spectral Density");
            testCase.assertEqual(size(spsd_no_rms.data, 2), expected_segments - 1);
        end

        function test_convert_to_bands(testCase)
            % Test converting SPSD to custom and standard bands
            m_spsd = testCase.spsd;
            
            % 1. Convert to millidecade
            mdec = convert_to_millidecade(m_spsd);
            testCase.assertTrue(isstruct(mdec));
            testCase.assertEqual(mdec.name, 'Millidecade Sound Pressure Spectral Density');
            testCase.assertEqual(size(mdec.data, 2), size(m_spsd.data, 2));
            testCase.verifyFalse(any(isnan(mdec.data(:))));

            % 2. Convert to decidecade
            ddec = convert_to_decidecade(m_spsd);
            testCase.assertTrue(isstruct(ddec));
            testCase.assertEqual(ddec.name, 'Decidecade Sound Pressure Spectral Density');
            testCase.assertEqual(size(ddec.data, 2), size(m_spsd.data, 2));
            testCase.verifyFalse(any(isnan(ddec.data(:))));

            % 3. Convert to third octave
            toct = convert_to_third_octave(m_spsd);
            testCase.assertTrue(isstruct(toct));
            testCase.assertEqual(toct.name, 'Third-Octave Sound Pressure Spectral Density');
            testCase.assertEqual(size(toct.data, 2), size(m_spsd.data, 2));
            testCase.verifyFalse(any(isnan(toct.data(:))));

            % 4. Convert to custom bands
            cust = convert_to_custom_bands(m_spsd, 20, 2, false);
            testCase.assertTrue(isstruct(cust));
            testCase.assertEqual(size(cust.data, 2), size(m_spsd.data, 2));
            testCase.verifyFalse(any(isnan(cust.data(:))));
        end

        function test_apply_calibration(testCase)
            time = 0:0.1:9.9;
            freq = linspace(10, 1000, numel(time));
            spsd_data = rand(numel(time), numel(freq));
            spsd.data = spsd_data;
            spsd.time = time;
            spsd.freq = freq;
            spsd.units = 'V^2/Hz';

            f_curve = rand(1, numel(freq))';
            fill_value = 0.0;

            calibrated_spsd = apply_calibration(spsd, [freq',f_curve], fill_value);

            testCase.assertTrue(isstruct(calibrated_spsd));
            testCase.assertEqual(calibrated_spsd.units, "Pa^2/Hz");
            testCase.assertEqual(size(calibrated_spsd.data), size(spsd.data));
            testCase.verifyLessThan(max(calibrated_spsd.data), max(spsd.data));
        end

        function test_freq_loss(testCase)
            fmin = minimum_frequency(20, 1500, 1700);
            testCase.assertEqual(fmin, 39.84375);
        end

        function test_spsdl(testCase)
            td_spsdl = sound_pressure_spectral_density_level(testCase.spsd);

            cc = datetime([ ...
                '2023-02-04 15:05:08.000'
                '2023-02-04 15:05:08.502'
                '2023-02-04 15:05:09.003'
                '2023-02-04 15:05:09.505'
                '2023-02-04 15:05:10.007'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
            cd_spsdl = [
                35.772061577247996	53.885951515794524	54.313180208582452	47.946115821970807	24.179883661640915;
                61.602936043528160	62.793717041169742	64.208043551995516	52.348715690865937	60.092868343630663;
                60.430110697601087	62.127736199381765	63.592162829344915	53.767504875422254	65.278154534206323;
                61.021732667297471	63.382134116643563	56.682598279488126	60.444088944449106	68.497868657610923;
                62.106024136858018	53.869424935555941	55.637077987595177	63.060370296217918	66.973991003254199;
            ];

            testCase.assertLessThan(abs(td_spsdl.data(1:5,1:5) - cd_spsdl), 1e-6);
            testCase.assertLessThan(abs(posixtime(td_spsdl.time(1:5))' - posixtime(cc)), 1e-3);
        end

        function test_averaging(testCase)
            td_spsdl = sound_pressure_spectral_density_level(testCase.spsd);

            octave = [3, 2];
            td_spsdl_mean = band_aggregate(td_spsdl, octave, 10, 100000);

            lbin = 30;
            td_spsdl_50 = time_aggregate(td_spsdl_mean, lbin, "median");

            cc = datetime([ ...
                '2023-02-04 15:05:23.000'
                '2023-02-04 15:05:53.000'
                '2023-02-04 15:06:23.000'
                '2023-02-04 15:06:53.000'
                '2023-02-04 15:07:23.000'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
            cd_spsdl_50 = [
                63.452642368099660	62.682385486447949	62.890011366673811	63.464187980216835	61.577799520029686;
                64.402822675991317	63.982645344031987	65.334381885187696	62.433832731393835	64.002193070627015;
                65.879232812149041	65.367523459295711	65.470768947676163	64.462081277290977	63.062659821297906;
                65.589916391755892	65.179564656250449	66.532335459481544	66.836901111597541	64.724432690177494;
                71.712416396688070	71.542148577084191	71.397701080122886	71.198748953358432	70.985008130759127
            ];

            testCase.assertLessThan(abs(td_spsdl_50.data(1:5,1:5) - cd_spsdl_50), 1e-5);
            testCase.assertLessThan(abs(posixtime(td_spsdl_50.time(1:5))' - posixtime(cc)), 1);
        end

        function test_fmax_warning(testCase)
            fn = 1000;
            fmax = 1500;
            adjusted_fmax = fmax_warning(fn, fmax);
            testCase.assertEqual(adjusted_fmax, fn);

            fmax = 500;
            adjusted_fmax = fmax_warning(fn, fmax);
            testCase.assertEqual(adjusted_fmax, fmax);

            testCase.verifyError(@() fmax_warning('not a number', fmax), 'MATLAB:validators:mustBeNumeric');
            testCase.verifyError(@() fmax_warning(fn, 'not a number'), 'MATLAB:validators:mustBeNumeric');
        end

        function test_validate_method(testCase)
            [method_name, method_arg] = validate_method("median");
            testCase.assertEqual(method_name, "median");
            testCase.assertEmpty(method_arg);

                [method_name, method_arg] = validate_method({"quantile", 0.25});
            testCase.assertEqual(method_name, "quantile");
            testCase.assertEqual(method_arg, 0.25);

            testCase.verifyError(@() validate_method('unsupported_method'), 'MHKiT:acoustics:validate_method');
            testCase.verifyError(@() validate_method(struct('unsupported_method', [])), 'MHKiT:acoustics:validate_method');
            testCase.verifyError(@() validate_method({'unsupported_method', 0.5}), 'MHKiT:acoustics:validate_method');
        end
    end
end