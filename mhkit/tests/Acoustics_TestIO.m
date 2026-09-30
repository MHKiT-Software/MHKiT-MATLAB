% filepath: c:\Github\MHKiT-Python\mhkit\tests\acoustics\test_io.m
classdef Acoustics_TestIO < matlab.unittest.TestCase
    properties
        datadir
        plotdir
    end

    methods (TestClassSetup)
        function setupOnce(testCase)
            testdir = fileparts(mfilename('fullpath'));
            testCase.plotdir = fullfile(testdir, 'plots');
            if ~isfolder(testCase.plotdir)
                mkdir(testCase.plotdir);
            end
            testCase.datadir = fullfile(testdir, '..', '..', 'examples', 'data', 'acoustics');
        end
    end

    methods (Test)

        function test_read_iclisten_metadata(testCase)
            file_name = fullfile(testCase.datadir, 'RBW_6661_20240601_053114.wav');
            metadata = read_iclisten(file_name);

            expected_metadata = struct(...
                'bits_per_sample', 24, ...
                'filename', 'RBW_6661_20240601_053114.wav', ...
                'peak_voltage', 3.0, ...
                'humidity', '24.0 % RH', ...
                'temperature', '8.6 deg C', ...
                'accelerometer', 'Acc(-980,-18,141)', ...
                'magnetometer', 'Mag(3603,3223,-598)', ...
                'count_at_peak_voltage', '8388608 = Max Count', ...
                'sequence_num', '2589798400000 = Seq #' ...
            );

            fields = fieldnames(expected_metadata);
            for i = 1:numel(fields)
                key = fields{i};
                expected_value = expected_metadata.(key);
                testCase.verifyTrue(isfield(metadata, key));
                if isnumeric(expected_value)
                    testCase.verifyEqual(metadata.(key), expected_value, 'AbsTol', 1e-6);
                else
                    testCase.verifyEqual(metadata.(key), expected_value);
                end
            end
        end

        function test_read_iclisten(testCase)
            file_name = fullfile(testCase.datadir, 'RBW_6661_20240601_053114.wav');
            td_orig = read_iclisten(file_name);
            td_wrap = read_hydrophone(file_name, 3, -177, 0, "2024-06-01T05:31:14");
            td_volt = read_iclisten(file_name, [], false);
            td_ovrrd = read_iclisten(file_name, -180, false);
            td_ovrrd2 = read_iclisten(file_name, -180, true);

            cc = datetime([...
                "2024-06-01T05:31:14.000000000"
                "2024-06-01T05:31:14.000001953"
                "2024-06-01T05:31:14.000003906"
                "2024-06-01T05:31:14.000005859"
                "2024-06-01T05:31:14.000007812"], 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS');

            cd_orig = [0.31546374, 0.30229832, 0.32229963, 0.3159701, 0.30356423];
            cd_volt = [0.0004456, 0.00042701, 0.00045526, 0.00044632, 0.0004288];
            cd_ovrrd = [0.44560438, 0.42700773, 0.45526033, 0.44631963, 0.42879587];

            testCase.verifyEqual(td_orig.data(1:5), cd_orig', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_orig.time(1:5)-cc', seconds(1e-9));

            testCase.verifyEqual(td_wrap.data(1:5), cd_orig', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_wrap.time(1:5)-cc', seconds(1e-9));

            testCase.verifyEqual(td_volt.data(1:5), cd_volt', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_volt.time(1:5)-cc', seconds(1e-9));

            testCase.verifyEqual(td_ovrrd.data(1:5), cd_ovrrd', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_ovrrd.time(1:5)-cc', seconds(1e-9));

            testCase.verifyEqual(td_ovrrd.data(1:5), td_ovrrd2.data(1:5), 'AbsTol', 1e-6);
        end

        function test_read_soundtrap(testCase)
            file_name = fullfile(testCase.datadir, '6247.230204150508.wav');
            td_orig = read_soundtrap(file_name, -177);
            td_wrap = read_hydrophone(file_name, 1, -177, 0, '2023-02-04T15:05:08');
            td_volt = read_soundtrap(file_name, []);

            cc = datetime([...
                "2023-02-04T15:05:08.000000000"
                "2023-02-04T15:05:08.000010416"
                "2023-02-04T15:05:08.000020832"
                "2023-02-04T15:05:08.000031249"
                "2023-02-04T15:05:08.000041665"], 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS');

            cd_orig = [0.929006, 0.929006, 0.929006, 0.929006, 1.01542517];
            cd_volt = [0.00131226, 0.00131226, 0.00131226, 0.00131226, 0.00143433];

            testCase.verifyEqual(td_orig.data(1:5), cd_orig', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_orig.time(1:5)-cc', seconds(1e-8));

            testCase.verifyEqual(td_wrap.data(1:5), cd_orig', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_wrap.time(1:5)-cc', seconds(1e-8));

            testCase.verifyEqual(td_volt.data(1:5), cd_volt', 'AbsTol', 1e-6);
            testCase.verifyLessThanOrEqual(td_volt.time(1:5)-cc', seconds(1e-8));
        end

        function test_calibration(testCase)
            file_name = fullfile(testCase.datadir, '6247.230204150508.wav');
            td_volt = read_soundtrap(file_name, []);
            td_spsd = sound_pressure_spectral_density(td_volt, td_volt.fs, 1);

            cal_name = fullfile(testCase.datadir, '6247_calibration.csv');
            calibration = readtable(cal_name);
            calibration = table2array(calibration);
            fill_Sf = calibration(end, 2); % use last value in curve for higher frequencies

            td_spsd_cal = apply_calibration(td_spsd, calibration, fill_Sf);

            cc = datetime([ ...
                '2023-02-04 15:05:08.000'
                '2023-02-04 15:05:08.502'
                '2023-02-04 15:05:09.003'
                '2023-02-04 15:05:09.505'
                '2023-02-04 15:05:10.007'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

            cd_spsd = [...
                0.000168347	0.010904234	0.012031451	0.002777228	1.16678E-05;
                0.004963596	0.006529418	0.009042888	0.000589351	0.00350582;
                0.000809366	0.001196485	0.001676293	0.000174536	0.002471435;
                0.000301447	0.000519099	0.000110993	0.000263904	0.00168587;
                0.000160878	2.41456E-05	3.62746E-05	0.000200416	0.000493508
                ];

            testCase.verifyEqual(td_spsd_cal.data(1:5,1:5), cd_spsd, 'AbsTol', 1e-6);
            testCase.verifyEqual(posixtime(td_spsd_cal.time(1:5)), posixtime(cc'), 'AbsTol', 1e-3);
        end

        function test_read_wispr_metadata(testCase)
            file_name = fullfile(testCase.datadir, 'WISPR_230825_003936.dat');
            metadata = read_wispr_metadata(file_name);

            expected_metadata = struct(...
                'version', 1.2, ...
                'time', "08:25:23:00:39:36", ...
                'instrument_id', "PERI_1", ...
                'location_id', "PWSPNE", ...
                'volts', 15.77, ...
                'blocks_free', 20.98, ...
                'file_size', 58575, ...
                'buffer_size', 16896, ...
                'samples_per_buffer', 8448, ...
                'sample_size', 2, ...
                'sampling_rate', 50000, ...
                'gain', 0, ...
                'decimation', 16, ...
                'adc_vref', 5, ...
                'file_length_sec', 299.904 ...
            );

            fields = fieldnames(expected_metadata);
            for i = 1:numel(fields)
                key = fields{i};
                expected_value = expected_metadata.(key);
                testCase.verifyTrue(isfield(metadata, key));
                if isnumeric(expected_value)
                    testCase.verifyEqual(metadata.(key), expected_value, 'AbsTol', 1e-4);
                else
                    testCase.verifyEqual(metadata.(key), expected_value);
                end
            end
        end

        function test_read_wispr(testCase)
            file_name = fullfile(testCase.datadir, 'WISPR_230825_003936.dat');
            out = read_wispr(file_name);

            testCase.verifyEqual(out.units, 'V');
            testCase.verifyEqual(out.fs, 50000);
            testCase.verifyEqual(out.peak_voltage, 5);
            testCase.verifyEqual(out.valid_min, -5);
            testCase.verifyEqual(out.valid_max, 5);

            expected_first_voltages = [-0.001678466796875, -0.001678466796875, -0.00152587890625, -0.0018310546875, -0.001068115234375]';
            testCase.verifyEqual(out.data(1:5), expected_first_voltages, 'AbsTol', 1e-9);

            cc = datetime([ ...
                "2023-08-25T00:39:36.000000"
                "2023-08-25T00:39:36.000020"
                "2023-08-25T00:39:36.000040"
                "2023-08-25T00:39:36.000060"
                "2023-08-25T00:39:36.000080"], 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

            testCase.verifyLessThanOrEqual(abs(out.time(1:5) - cc), seconds(1e-6));
        end

        function test_export_audio(testCase)
            file_name = fullfile(testCase.datadir, '6247.230204150508.wav');
            td_pressure = read_soundtrap(file_name, -177);

            % Create output filename in plots directory
            output_filename = fullfile(testCase.plotdir, 'test_export_audio');

            % Test export_audio with default gain
            export_audio(output_filename, td_pressure);

            % Verify the output file was created
            expected_wav_file = strcat(output_filename, '.wav');
            testCase.verifyTrue(exist(expected_wav_file, 'file') == 2, ...
                'Export audio should create a WAV file');

            % Verify the audio file can be read back
            [audio_data, fs_read] = audioread(expected_wav_file);
            testCase.verifyEqual(fs_read, td_pressure.fs, ...
                'Exported audio sample rate should match input');
            testCase.verifyEqual(length(audio_data), length(td_pressure.data), ...
                'Exported audio length should match input');

            % Clean up test file
            if exist(expected_wav_file, 'file')
                delete(expected_wav_file);
            end

            % Test export_audio with resampling multiplier
            resampled_output_filename = fullfile(testCase.plotdir, 'test_export_audio_resampled');
            export_audio(resampled_output_filename, td_pressure, [], 1, 2);
            expected_resampled_wav = strcat(resampled_output_filename, '.wav');
            testCase.verifyTrue(exist(expected_resampled_wav, 'file') == 2, ...
                'Export resampled audio should create a WAV file');
            [audio_data_res, fs_res] = audioread(expected_resampled_wav);
            testCase.verifyEqual(fs_res, td_pressure.fs);
            testCase.verifyEqual(length(audio_data_res), floor(length(td_pressure.data)/2));
            if exist(expected_resampled_wav, 'file')
                delete(expected_resampled_wav);
            end
        end
    end
end
