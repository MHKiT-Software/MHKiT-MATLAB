classdef Wave_TestIO < matlab.unittest.TestCase
    
    methods (Test)
        
        % Realtime data
        function test_read_NDBC_realtime_met(testCase)
            Obj.expected_columns_metRT = struct('WDIR',{},'units',{},'WSPD',{},'GST',{},'WVHT',{},'DPD',{},'APD',{},'MWD',{},'PRES',{},'ATMP',{},'WTMP',{},'DEWP',{},'VIS',{},'PTDY',{},'TIDE',{},'time',{});
            Obj.expected_units_metRT = struct('WDIR',{"degT"},'WSPD',{"m/s"},'GST',{"m/s"},'WVHT',{"m"},'DPD',{"sec"},'APD',{"sec"},'MWD',{"degT"},'PRES',{"hPa"},'ATMP',{"degC"},'WTMP',{"degC"},'DEWP',{"degC"},'VIS',{"nmi"},'PTDY',{"hPa"},'TIDE',{"ft"});

            data = read_NDBC_file("../../examples/data/wave/46097.txt");
            datarm1 = rmfield(data,'units');
            datarm2 = rmfield(datarm1,'time');
            expected_index0 = datestr(datenum(2019,4,2,13,50,0));
            
            assertEqual(testCase,fieldnames(data),fieldnames(Obj.expected_columns_metRT));      
            assertEqual(testCase,datestr(data.time(1)), expected_index0);
            assertEqual(testCase,size(getfield(datarm2,'WDIR')), [6490 1]);
            assertEqual(testCase,size(getfield(datarm2,'TIDE')), [6490 1]);
            assertEqual(testCase,data.units,Obj.expected_units_metRT);
        end

        % Historical data
        function test_read_NDBC_historical_met(testCase)
            Obj.expected_columns_metH = struct('WDIR',{},'units',{},'WSPD',{},'GST',{},'WVHT',{},'DPD',{},'APD',{},'MWD',{},'PRES',{},'ATMP',{},'WTMP',{},'DEWP',{},'VIS',{},'TIDE',{},'time',{});
            Obj.expected_units_metH = struct('WDIR',{"degT"},'WSPD',{"m/s"},'GST',{"m/s"},'WVHT',{"m"},'DPD',{"sec"},'APD',{"sec"},'MWD',{"deg"},'PRES',{"hPa"},'ATMP',{"degC"},'WTMP',{"degC"},'DEWP',{"degC"},'VIS',{"nmi"},'TIDE',{"ft"});
            
            data = read_NDBC_file('../../examples/data/wave/46097h201908qc.txt');
            datarm1 = rmfield(data,'units');
            datarm2 = rmfield(datarm1,'time');
            expected_index0 = datestr(datenum(2019,8,1,0,0,0));
            
            assertEqual(testCase,fieldnames(data),fieldnames(Obj.expected_columns_metH));
            assertEqual(testCase,datestr(data.time(1)), expected_index0);
            assertEqual(testCase,size(getfield(datarm2,'WDIR')), [4464 1]);
            assertEqual(testCase,size(getfield(datarm2,'TIDE')), [4464 1]);
            assertEqual(testCase,data.units,Obj.expected_units_metH);
        end

        % Spectral data
        function test_read_NDBC_spectral(testCase)
            data = read_NDBC_file("../../examples/data/wave/data.txt");
            assertEqual(testCase,size(data.spectra),[743 47]);
        end
    end
end