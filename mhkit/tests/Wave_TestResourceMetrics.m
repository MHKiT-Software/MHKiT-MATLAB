classdef Wave_TestResourceMetrics < matlab.unittest.TestCase

    methods (Test)

        function test_kfromw(testCase)
            relative_file_name = '../../examples/data/wave/ValData1.json'; % filename in JSON extension
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string

            Data1_valdata1_table = struct2table(valdata1.Data1,'AsArray',true);
            Data1_valdata1 = table2struct(Data1_valdata1_table,'ToScalar',true);
            Data11_valdata1_table = struct2table(valdata1.Data11,'AsArray',true);
            Data11_valdata1 = table2struct(Data11_valdata1_table,'ToScalar',true);

            f1 = valdata1.Data1.w/(2*pi);
            f2 = valdata1.Data11.w/(2*pi);
            f1 = f1';
            f2 = f2';
            h1 = Data1_valdata1.h;
            h2 = Data11_valdata1.h;
            rho1 = Data1_valdata1.rho;
            rho2 = Data11_valdata1.rho;
            expected1 = valdata1.Data1.k;
            expected2 = valdata1.Data11.k;
            calculated1 = wave_number(f1, h1, "rho", rho1);
            calculated2 = wave_number(f2, h2, "rho",rho2);
            error1 = sum(((expected1-calculated1.frequency)^2)); % SSE
            error2 = sum(((expected2-calculated2.frequency)^2)); % SSE
            assertLessThan(testCase,error1, 1e-6);
            assertLessThan(testCase,error2, 1e-6);
        end

        function test_moments(testCase)
            format long
            relative_file_name = '../../examples/data/wave/ValData2.mat';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = load(full_file_name);
            Valdata = data.CalcSpecCheckData;

            H5SP = Valdata.H5sP;
            H10SP = Valdata.H10sP;
            AH1 = Valdata.AH1;
            AH6 = Valdata.AH6;
            CDIP1 = Valdata.CDiP1;
            CDIP6 = Valdata.CDiP6;

            S1 = struct('spectrum',H5SP.PowSpecHanWin','type','TimeSeries','frequency',H5SP.freq);
            S2 = struct('spectrum',H10SP.PowSpecHanWin','type','TimeSeries','frequency',H10SP.freq);
            S3 = struct('spectrum',AH1.PowSpecHanWin,'type','TimeSeries','frequency',AH1.freq);
            S4 = struct('spectrum',AH6.PowSpecHanWin,'type','TimeSeries','frequency',AH6.freq);
            S5 = struct('spectrum',CDIP1.PowSpec,'type','TimeSeries','frequency',CDIP1.freq);
            S6 = struct('spectrum',CDIP6.PowSpec,'type','TimeSeries','frequency',CDIP6.freq);

            expected1 = H5SP.Moment;
            expected2 = H10SP.Moment;
            expected3 = AH1.Moment;
            expected4 = AH6.Moment;
            expected5 = CDIP1.Moment;
            expected6 = CDIP6.Moment;

            for i = 0:5
                calculated1 = frequency_moment(S1, i-1);
                assertEqual(testCase,expected1(i+1),calculated1,'AbsTol',0.01);
                calculated2 = frequency_moment(S2, i-1);
                assertEqual(testCase,expected2(i+1),calculated2,'AbsTol',0.01);
                calculated3 = frequency_moment(S3, i-1);
                assertEqual(testCase,expected3(i+1),calculated3,'AbsTol',0.01);
                calculated4 = frequency_moment(S4, i-1);
                assertEqual(testCase,expected4(i+1),calculated4,'AbsTol',0.01);
                calculated5 = frequency_moment(S5, i-1,CDIP1.freqBinWidth);
                assertEqual(testCase,expected5(i+1),calculated5,'AbsTol',0.01);
                calculated6 = frequency_moment(S6, i-1,CDIP6.freqBinWidth);
                assertEqual(testCase,expected6(i+1),calculated6,'AbsTol',0.01);
            end
        end

       function test_metrics_HsP(testCase)

            relative_file_name = '../../examples/data/wave/ValData2.mat';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = load(full_file_name);
            Valdata = data.CalcSpecCheckData;
            H5SP = Valdata.H5sP;
            S1 = struct('spectrum',H5SP.PowSpecHanWin','type','TimeSeries','frequency',H5SP.freq);

            % Hm0
            expected = H5SP.waveMoments.Hm0;
            calculated = significant_wave_height(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Te
            expected = H5SP.waveMoments.Te;
            calculated = energy_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % T0
            expected = H5SP.waveMoments.T0;
            calculated = average_zero_crossing_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tc
            expected = H5SP.waveMoments.Tc;
            calculated = average_crest_period(S1)^2;
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tm
            expected = sqrt(H5SP.waveMoments.Tm);
            calculated = average_wave_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tp
            expected = H5SP.waveMoments.Tp;
            calculated = peak_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % e
            expected = H5SP.waveMoments.e;
            calculated = spectral_bandwidth(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % v
            expected = H5SP.waveMoments.v;
            calculated = spectral_width(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);
       end

       function test_metrics_AH(testCase)

            relative_file_name = '../../examples/data/wave/ValData2.mat';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = load(full_file_name);
            Valdata = data.CalcSpecCheckData;
            AH1 = Valdata.AH1;
            S1 = struct('spectrum',AH1.PowSpecHanWin,'type','TimeSeries','frequency',AH1.freq);

            % Hm0
            expected = AH1.waveMoments.Hm0;
            calculated = significant_wave_height(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Te
            expected = AH1.waveMoments.Te;
            calculated = energy_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % T0
            expected = AH1.waveMoments.T0;
            calculated = average_zero_crossing_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tc
            expected = AH1.waveMoments.Tc;
            calculated = average_crest_period(S1)^2;
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tm
            expected = sqrt(AH1.waveMoments.Tm);
            calculated = average_wave_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tp
            expected = AH1.waveMoments.Tp;
            calculated = peak_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % e
            expected = AH1.waveMoments.e;
            calculated = spectral_bandwidth(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % v
            expected = AH1.waveMoments.v;
            calculated = spectral_width(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);
       end

       function test_metrics_CDIP1(testCase)

            relative_file_name = '../../examples/data/wave/ValData2.mat';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = load(full_file_name);
            Valdata = data.CalcSpecCheckData;
            CDIP1 = Valdata.CDiP1;
            S1 = struct('spectrum',CDIP1.PowSpec,'type','TimeSeries','frequency',CDIP1.freq);

            % Hm0
            expected = CDIP1.waveMoments.Hm0;
            calculated = significant_wave_height(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Te
            expected = CDIP1.waveMoments.Te;
            calculated = energy_period(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % T0
            expected = CDIP1.waveMoments.T0;
            calculated = average_zero_crossing_period(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tc
            expected = CDIP1.waveMoments.Tc;
            calculated = average_crest_period(S1,CDIP1.freqBinWidth)^2;
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tm
            expected = sqrt(CDIP1.waveMoments.Tm);
            calculated = average_wave_period(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);

            % Tp
            expected = CDIP1.waveMoments.Tp;
            calculated = peak_period(S1);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % e
            expected = CDIP1.waveMoments.e;
            calculated = spectral_bandwidth(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.001);

            % v
            expected = CDIP1.waveMoments.v;
            calculated = spectral_width(S1,CDIP1.freqBinWidth);
            error = abs(expected-calculated)/expected;
            assertLessThan(testCase,error, 0.01);
        end

        function test_plot_elevation_timeseries(testCase)

            relative_file_name = '../../examples/data/wave/ValData2.mat';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = load(full_file_name);
            Valdata = data.CalcSpecCheckData;
            H5sP = Valdata.H5sP;

            df = 0.01/(2*pi);
            Trep = 1/df;
            time = 0:0.062838:Trep;

            elevation = H5sP.TimeSeries;
            sample_rate = H5sP.Fs;
            NFFT = H5sP.NFFT;
            H5sP.S = elevation_spectrum(elevation',sample_rate,NFFT,time);

            filename = 'wave_plot_elevation_timeseries.png';
            if isfile(filename)
                delete(filename);
            end

            wave_elevation = struct('time',time','elevation',H5sP.S.spectrum);

            df = 0.01/(2*pi);
            Trep = 1/df;
            time = 0:0.12566:Trep;
            wave_elevation.time = time';

            plot_elevation_timeseries(wave_elevation,"savepath",filename);
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_environmental_contour(testCase)

            % assumeFail(testCase, "Not compatible with latest MHKIT-Python")

            relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);

            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Te_table = struct2table(valdata1.Te,'AsArray',true);
            Te = table2array(Te_table);
            Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
            Hm0 = table2array(Hm0_table);

            filter = Hm0 < 20;
            Hm0 = Hm0(filter);
            Te = Te(filter);
            [row, col] = find(~isnan(Te));
            Hm0 = Hm0(col);
            Te = Te(col);
            [row, col] = find(~isnan(Hm0));
            Hm0 = Hm0(col);
            Te = Te(col);

            time_str = Hm0_table.Properties.VariableNames;

            time1 = str2num(erase(time_str{1},'x'));
            time2 = str2num(erase(time_str{2},'x'));

            dt = (time2-time1)/1000.;
            time_R = 100;

            contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');

            relative_file_name= '../../examples/data/wave/Hm0_Te_contours_46022.csv';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            expected_contours = readmatrix(full_file_name);

            Hm0_expected = expected_contours(:,1);
            Te_expected = expected_contours(:,2);

            assertEqual(testCase,contour.contour1,Hm0_expected','RelTol',0.01);
            assertEqual(testCase,contour.contour2,Te_expected','RelTol',0.01);


        end

        function test_plot_environmental_contour(testCase)

            %assumeFail(testCase, "Not compatible with latest MHKIT-Python")

            relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);

            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Te_table = struct2table(valdata1.Te,'AsArray',true);
            Te = table2array(Te_table);
            Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
            Hm0 = table2array(Hm0_table);

            filter = Hm0 < 20;
            Hm0 = Hm0(filter);
            Te = Te(filter);
            [row, col] = find(~isnan(Te));
            Hm0 = Hm0(col);
            Te = Te(col);
            [row, col] = find(~isnan(Hm0));
            Hm0 = Hm0(col);
            Te = Te(col);

            time_str = Hm0_table.Properties.VariableNames;

            time1 = str2num(erase(time_str{1},'x'));
            time2 = str2num(erase(time_str{2},'x'));

            dt = (time2-time1)/1000.;
            time_R = 100;

            contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');

            filename = 'wave_plot_env_contour.png';
            if isfile(filename)
                delete(filename);
            end


            plot_environmental_contours(Te, Hm0,contour.contour2,contour.contour1,"savepath",filename...
                ,"x_label",...
                'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
                "contour_label",'100 Year Contour');
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_plot_environmental_contour_multiyear(testCase)

            assumeFail(testCase, "Not compatible with latest MHKIT-Python")
            % not sure about why this test exists...return period has to be float or
            % int and cannot be a list...
            relative_file_name= '../../examples/data/wave/Hm0_Te_46022.json';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);

            fid = fopen(full_file_name); % Opening the file
            raw = fread(fid,inf); % Reading the contents
            str = char(raw'); % Transformation
            fclose(fid); % Closing the file
            valdata1 = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Te_table = struct2table(valdata1.Te,'AsArray',true);
            Te = table2array(Te_table);
            Hm0_table = struct2table(valdata1.Hm0,'AsArray',true);
            Hm0 = table2array(Hm0_table);

            filter = Hm0 < 20;
            Hm0 = Hm0(filter);
            Te = Te(filter);
            [row, col] = find(~isnan(Te));
            Hm0 = Hm0(col);
            Te = Te(col);
            [row, col] = find(~isnan(Hm0));
            Hm0 = Hm0(col);
            Te = Te(col);

            time_str = Hm0_table.Properties.VariableNames;

            time1 = str2num(erase(time_str{1},'x'));
            time2 = str2num(erase(time_str{2},'x'));

            dt = (time2-time1)/1000.;
            time_R = [100, 120, 130];

            contour = environmental_contours(Hm0, Te, dt, time_R, 'PCA');

            filename = 'wave_plot_env_contour_multiyear.png';
            if isfile(filename)
                delete(filename);
            end


            plot_environmental_contours(Te, Hm0,contour.contour2,contour.contour1,"savepath",filename...
                ,"x_label",...
                'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
                "contour_label",{'100 Year Contour','120 Year Contour','130 Year Contour'});
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_chakrabarti(testCase)
            D = 5;
            H = 8;
            lambda_w = 200;
            filename = 'chakrabarti.png';
            if isfile(filename)
                delete(filename);
            end
            plot_chakrabarti(H,lambda_w,D,'savepath',filename)
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end

        function test_wave_length(testCase)
            k=[1,2,10,3];
            l_expected = (2.*3.14)./k;


            l_calculated = wave_length(k);
            assertEqual(testCase,l_expected,l_calculated,'RelTol',0.01);
        end


         function test_depth_regime(testCase)
             expected = [1,1,0,1];
             l_vector=[1,2,10,3];

             h = 10;

             calculated = depth_regime(l_vector,h);
             assertEqual(testCase,expected,calculated);
         end

        function test_wave_celerity(testCase)
            % small change in f will give similar value cg
            f=linspace(20.0001,20.0005,5);

            % Choose index to spike at. cg spike is inversly proportional to k
            k_tmp=[1, 1, 0.5, 1, 1];
            k.values = k_tmp;
            k.frequency=f;

            % all shallow
            cg_shallow1 = wave_celerity(k,0.0001,"depth_check",py.True);
            cg_shallow2 = wave_celerity(k, 0.0001,"depth_check",py.False);
            assertEqual(testCase,cg_shallow1, cg_shallow2);

            x = (3.14.*f)./k.values;
            % all deep
            cg = wave_celerity(k, 1000,"depth_check",py.True);
            assertEqual(testCase,x,cg.values,'RelTol',0.01);
        end

        function test_energy_flux_deep(testCase)
            omega = [0.1:0.01:3.5];
            f = omega./(2*pi);
            Hs = 2.5;
            Tp = 8 ;
            % Dependent on mhkit.resource.BS spectrum
            S = jonswap_spectrum(f,Tp,Hs);
            Te = energy_period(S);
            Hm0 = significant_wave_height(S);
            rho=1025;
            g=9.80665;

            coeff = rho*(g^2)/(64*pi);
            J = coeff*(Hm0^2)*Te;

            h=-1; % not used when deep=True
            J_calc = energy_flux(S, h, "deep",py.True);

            assertEqual(testCase,J_calc,J,'RelTol',0.01);
        end
    end
end
