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
            
            assertEqual(testCase,size(data.spectrum),[47 743]);
        end
        
        function test_ndbc_available_data(testCase)
            data = NDBC_available_data('swden', 'buoy_number','46029');
            columns = data.Properties.VariableNames;
            columns_expected = [{'Station_id'},{'year'},{'file'}];
            assertEqual(testCase,columns,columns_expected);
            number_years = length(data.year);
            expected_years = [1996:1996+(number_years-1)]';
            assertEqual(testCase,str2num(char(data.year)),expected_years);
            
            h = height(data);
            w = width(data);
            assertEqual(testCase,[h,w],[number_years,3]);
            
         
        end
        
        function test_ndbc_request_data(testCase)
            filenames=["46042w1996.txt.gz";... 
                        "46029w1997.txt.gz";.... 
                        "46029w1998.txt.gz"];
            ndbc_data = NDBC_request_data('swden', filenames);
            ndbc_data = struct2table(ndbc_data.ID_46042.year_1996,'AsArray',true);
            file = gunzip('../../examples/data/wave/46042w1996.txt.gz');
            expected_data = readmatrix(file{1});
            temp = table2array(ndbc_data(:,6));
            data_array = table2array(ndbc_data(:,[1,2,3,4]));
            temp = cellfun(@transpose,temp,'UniformOutput',false);
            cat_data = cat(2,data_array{:},temp{:});
            
            assertEqual(testCase,cat_data,expected_data(2:end,:));
            
         
        end
    end
end