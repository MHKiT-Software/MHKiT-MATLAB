classdef Loads_TestExtreme < matlab.unittest.TestCase

    methods (Test, TestTags = {'DebuggingActions'})

        function test_mhkit_import(testCase)
            py.importlib.import_module('mhkit');
            assertTrue(testCase, true);
        end

    end

    methods (Test)

        function test_mler_coefficients(testCase)
            % create inputs and load validation data
            fpath = '../../examples/data/loads/mler.csv';
            validation = readtable(fpath);
            wave_freq = linspace(0,1,500);
            js = jonswap_spectrum(wave_freq,15.1,9);
            response_desired = 1;
            RAO = validation.RAO;
            % execute function
            mler = mler_coefficients(RAO, js, response_desired);
            % assertions
            assertEqual(testCase, mler.conditioned_spectrum, validation.Res_Spec, 'RelTol',0.005)
            assertEqual(testCase, mler.phase, validation.phase, 'RelTol',0.001)
        end

        function test_mler_simulation(testCase)
            % validation data
            T = linspace(-150,150,301);
            X = linspace(-300,300,601);
            sim = mler_simulation();
            
            assertEqual(testCase, sim.T, T)
            assertEqual(testCase, sim.X, X)
        end
        
        function test_mler_wave_amp_normalize(testCase)
            fpath = '../../examples/data/loads/mler.csv';
            validation = readtable(fpath);
            mler.conditioned_spectrum = validation.Res_Spec;
            mler.phase = validation.phase;
            wave_freq = linspace(0,1,500);
            mler.frequency = wave_freq';
            
            k = wave_number(wave_freq, 70);
            k.values = fillmissing(k.values,'constant',0);
            sim = mler_simulation();
            mler_norm = mler_wave_amp_normalize(4.5*1.9, mler, sim, k.values);
            
            assertEqual(testCase, mler_norm.conditioned_spectrum, validation.Norm_Spec, 'AbsTol', 0.002)
        end
        
        function test_mler_export_time_series(testCase)
            fpath = '../../examples/data/loads/mler_ts.csv';
            validation = readtable(fpath);
            fpath2 = '../../examples/data/loads/mler.csv';
            normed = readtable(fpath2);
            mler.conditioned_spectrum = normed.Norm_Spec;
            mler.phase = normed.phase;
            wave_freq = linspace(0,1,500);
            mler.frequency = wave_freq';
            RAO = normed.RAO;
            k = wave_number(wave_freq, 70);
            k.values = fillmissing(k.values,'constant',0); 
            sim = mler_simulation();
            mler_ts = mler_export_time_series(RAO, mler, sim, k.values);
            
            assertEqual(testCase, mler_ts.linear_response, validation.LinearResponse, 'AbsTol', 0.00005)
        end
        
        function test_global_peaks(testCase)
            fpath = '../../examples/data/loads/sine_wave.csv';
            validation = readtable(fpath);
            Fs = 200;                   % samples per second
            dt = 1/Fs;                   % seconds per sample
            StopTime = 0.25;             % seconds
            t = (0:dt:StopTime-dt)';     % seconds
            Fc = 60;                     % hertz
            x = cos(2*pi*Fc*t);
            [t_peaks, peaks] = global_peaks(t, x);

            assertEqual(testCase, t_peaks, validation.t_peaks', 'AbsTol', 0.00005)
            assertEqual(testCase, peaks, validation.peaks', 'AbsTol', 0.00005)
        end
        
        function test_number_of_short_term_peaks(testCase)
            n = 10;
            t = 100;
            t_st = 50;
            n_st = number_of_short_term_peaks(n, t, t_st);

            assertEqual(testCase, n_st, 5)
        end

        function test_block_maxima(testCase)
            load('../../examples/data/loads/block_maxima.mat', 'block_max');
            load('../../examples/data/loads/example_qoi.mat', 'qoi')
            % time in seconds
            time = linspace(0, 10800, 2*10801+1);
            % split time series into 10-minute (600s) segments
            bm = block_maxima(time, qoi', 600);

            assertEqual(testCase, bm, block_max)
        end

        function test_ste_block_maxima_gev(testCase)
            load('../../examples/data/loads/block_maxima.mat', 'block_max');
            load('../../examples/data/loads/ste.mat', 'ste_gev_pdf')
            load('../../examples/data/loads/ste.mat', 'ste_gev_cdf')
            x = linspace(0,3,1000);
            stegev_pdf = ste_block_maxima_gev(block_max, x, "pdf");
            stegev_cdf = ste_block_maxima_gev(block_max, x, "cdf");

            assertEqual(testCase, stegev_pdf, ste_gev_pdf)
            assertEqual(testCase, stegev_cdf, ste_gev_cdf)
        end
        
        
    end
end