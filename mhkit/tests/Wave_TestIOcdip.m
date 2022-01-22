classdef Wave_TestIOcdip < matlab.unittest.TestCase
    %WAVE_TESTIOCDIP Testing CDIP data querying and formatting
    %   The Coastal Data Information Program (CDIP) measures, analyzes,
    %   archives and disseminates coastal environment data.
    %   https://cdip.ucsd.edu/
    
%     properties
%         Property1
%     end
    
    methods (Test)
%         function obj = Wave_TestIOcdip(inputArg1,inputArg2)
%             %WAVE_TESTIOCDIP Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
%         
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end

        function test_request_parse_workflow_multiyear(testCase)
            station_number = '067';
            year1 = 2011;
            year2 = 2013;
            years = [year1, year2];
            parameters = {'waveHs', 'waveMeanDirection', 'waveA1Value'};
            data = cdip_request_parse_workflow( ...
                'station_number', station_number, ...
                'years', years, ...
                'parameters', parameters);

            assertEqual(testCase,1,1);

%             expected_index0 = datetime(year1, 1, 1);
%             expected_index_final = datetime(year2, 12, 30); % last data on 30th
%
%             wave1D = data['data']['wave'];
%             assertEqual(testCase,wave1D.index[0].floor('d').to_pydatetime(), expected_index0);
%             assertEqual(testCase,wave1D.index[-1].floor('d').to_pydatetime(), expected_index_final);
% 
%             for key,wave2D  in data['data']['wave2D'].items():
%                 assertEqual(wave2D.index[0].floor('d').to_pydatetime(), expected_index0);
%                 assertEqual(wave2D.index[-1].floor('d').to_pydatetime(), expected_index_final);
%             end
        end
    end
end

