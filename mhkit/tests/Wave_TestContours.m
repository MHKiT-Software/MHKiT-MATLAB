classdef Wave_TestContours < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function test_samples_contour(testCase)
            file_loc = "../../examples/data/wave/WDRT_caluculated_countours.json";
            str = fileread(file_loc); % dedicated for reading files as text 
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string 
            te_samples = [10, 15, 20];
            hs_samples_0 = [8.56637939, 9.27612515, 8.70427774];
            hs_contour = data.gaussian_x1;
            te_contour = data.gaussian_x2;
            hs_samples = samples_contour(te_samples, te_contour, hs_contour);
            
            assertEqual(testCase, hs_samples, hs_samples_0, 'AbsTol', 0.0005)
        end

        function test_samples_full_seastate(testCase)
            hs_0 = [5.91760129, 4.55185088, 1.41144991, 12.64443154, 7.89753791, 0.93890797];
            te_0 = [14.24199604, 8.25383556, 6.03901866, 16.9836369, 9.51967777, 3.46969355];
            w_0 = [2.18127398e-01,2.18127398e-01,2.18127398e-01,2.45437862e-07,2.45437862e-07,2.45437862e-07];

            file_loc = "../../examples/data/wave/Hm0_Te_46022.json";
            str = fileread(file_loc); % dedicated for reading files as text 
            data = jsondecode(str); % Using the jsondecode function to parse JSON from string
            Hm0 = cell2mat(struct2cell(data.Hm0));
            qc = find(Hm0 < 20);
            Te = cell2mat(struct2cell(data.Te));
            dt_ss = 3600;
            points_per_interval = 3;
            return_periods = [50, 100];
            py.numpy.random.seed(int8(0));
            [hs, te, w] = samples_full_seastate(Hm0(qc), Te(qc), points_per_interval, return_periods, dt_ss, "PCA", 250);
            assertEqual(testCase, hs, hs_0, 'AbsTol',0.0005)
            assertEqual(testCase, te, te_0, 'AbsTol',0.0005)
            assertEqual(testCase, w, w_0, 'AbsTol',0.0005)
        end
    end
    
end