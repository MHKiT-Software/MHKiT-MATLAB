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
    end
end  




        