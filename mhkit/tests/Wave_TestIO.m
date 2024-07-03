classdef Wave_TestIO < matlab.unittest.TestCase

    methods (Test)

        % Realtime data
        function test_read_NDBC_realtime_met(testCase)
            Obj.expected_columns_metRT = struct('WDIR',{},'units',{},'WSPD',{},'GST',{},'WVHT',{},'DPD',{},'APD',{},'MWD',{},'PRES',{},'ATMP',{},'WTMP',{},'DEWP',{},'VIS',{},'PTDY',{},'TIDE',{},'time',{});
            Obj.expected_units_metRT = struct('WDIR',{"degT"},'WSPD',{"m/s"},'GST',{"m/s"},'WVHT',{"m"},'DPD',{"sec"},'APD',{"sec"},'MWD',{"degT"},'PRES',{"hPa"},'ATMP',{"degC"},'WTMP',{"degC"},'DEWP',{"degC"},'VIS',{"nmi"},'PTDY',{"hPa"},'TIDE',{"ft"});

            relative_file_name = "../../examples/data/wave/46097.txt";
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_NDBC_file(full_file_name);
            datarm1 = rmfield(data,'units');
            datarm2 = rmfield(datarm1,'time');
            expected_index0 = posixtime(datetime(2019,4,2,13,50,0));

            assertEqual(testCase,fieldnames(data),fieldnames(Obj.expected_columns_metRT));
            assertEqual(testCase,data.time(1),expected_index0);
            assertEqual(testCase,size(getfield(datarm2,'WDIR')), [6490 1]);
            assertEqual(testCase,size(getfield(datarm2,'TIDE')), [6490 1]);
            assertEqual(testCase,data.units,Obj.expected_units_metRT);
        end

        % Historical data
        function test_read_NDBC_historical_met(testCase)
            Obj.expected_columns_metH = struct('WDIR',{},'units',{},'WSPD',{},'GST',{},'WVHT',{},'DPD',{},'APD',{},'MWD',{},'PRES',{},'ATMP',{},'WTMP',{},'DEWP',{},'VIS',{},'TIDE',{},'time',{});
            Obj.expected_units_metH = struct('WDIR',{"degT"},'WSPD',{"m/s"},'GST',{"m/s"},'WVHT',{"m"},'DPD',{"sec"},'APD',{"sec"},'MWD',{"deg"},'PRES',{"hPa"},'ATMP',{"degC"},'WTMP',{"degC"},'DEWP',{"degC"},'VIS',{"nmi"},'TIDE',{"ft"});

            relative_file_name = '../../examples/data/wave/46097h201908qc.txt';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_NDBC_file(full_file_name);
            datarm1 = rmfield(data,'units');
            datarm2 = rmfield(datarm1,'time');
            expected_index0 = posixtime(datetime(2019,8,1,0,0,0));

            assertEqual(testCase,fieldnames(data),fieldnames(Obj.expected_columns_metH));
            assertEqual(testCase,data.time(1),expected_index0);
            assertEqual(testCase,size(getfield(datarm2,'WDIR')), [4464 1]);
            assertEqual(testCase,size(getfield(datarm2,'TIDE')), [4464 1]);
            assertEqual(testCase,data.units,Obj.expected_units_metH);
        end

        % Spectral data
        function test_read_NDBC_spectral(testCase)
            relative_file_name = "../../examples/data/wave/data.txt";
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            data = read_NDBC_file(full_file_name);

            assertEqual(testCase,size(data.spectrum),[47 743]);
        end

        function test_ndbc_available_data(testCase)
            data = NDBC_available_data('swden', 'buoy_number','46029');
            columns = fieldnames(data);
            columns_expected = [{'Station_id'};{'year'};{'file'}];
            assertEqual(testCase,columns,columns_expected);
            unique_years = unique(data.year);
            number_years = length(unique_years);
            expected_years = [1996:1996+(number_years-1)]';
            assertEqual(testCase,unique_years,expected_years);
        end

        function test_ndbc_request_data(testCase)
            filenames=["46042w1996.txt.gz";...
                        "46029w1997.txt.gz";....
                        "46029w1998.txt.gz"];
            ndbc_data = NDBC_request_data('swden', filenames);
            ndbc_data = struct2table(ndbc_data.ID_46042.year_1996,'AsArray',true);
            relative_file_name = '../../examples/data/wave/46042w1996.txt.gz';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            file = gunzip(full_file_name);
            expected_data = readmatrix(file{1});
            temp = table2array(ndbc_data(:,7));
            data_array = table2array(ndbc_data(:,[1,2,3,4]));
            temp = cellfun(@transpose,temp,'UniformOutput',false);
            cat_data = cat(2,[data_array{:}],temp{:});

            assertEqual(testCase,cat_data,expected_data(2:end,:));


        end

        function test_swan_read_table(testCase)
            relative_file_name = '../../examples/data/wave/SWAN/SWANOUT.DAT';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            delimiterIn = ' ';
            mystructure = importdata(full_file_name,delimiterIn);
            vars = string(strsplit(mystructure.textdata{5},' '));
            vars = vars(2:end-1);
            expected = table2struct(array2table(mystructure.data,'VariableNames',vars),'ToScalar',true);
            data = swan_read_table(full_file_name);

            assertEqual(testCase,data.Hsig,expected.Hsig);
        end

        function test_swan_read_block(testCase)
            relative_file_name = '../../examples/data/wave/SWAN/SWANOUT.DAT';
            relative_file_name2 = '../../examples/data/wave/SWAN/SWANOUTBlock.DAT';
            full_file_name = fullfile(fileparts(mfilename('fullpath')), relative_file_name);
            full_file_name2 = fullfile(fileparts(mfilename('fullpath')), relative_file_name2);
            delimiterIn = ' ';
            mystructure = importdata(full_file_name,delimiterIn);
            vars = string(strsplit(mystructure.textdata{5},' '));
            vars = vars(2:end-1);
            expected = table2struct(array2table(mystructure.data,'VariableNames',vars),'ToScalar',true);
            data = swan_read_block(full_file_name2);

            assertEqual(testCase,sum(sum(data.Significant_wave_height.values)),sum(expected.Hsig),'RelTol',0.001);
        end

        % WPTO multiple locations
        function test_WPTO_point_multiloc(testCase)

            assumeFail(testCase, "API key usage saturated - temporarily disabling")
            % Error using matlab.internal.webservices.HTTPConnector/copyContentToByteArray (line 373)
            % The server returned the status 503 with message "Service Unavailable" in response to the request to URL
            % https://developer.nrel.gov/api/hsds/?api_key=3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf&domain=%2Fnrel%2FUS_wave%2Fvirtual_buoy%2FWest_Coast%2FWest_Coast_virtual_buoy_2010.h5.
            api_key = '3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf';
            hindcast_data = request_wpto('1-hour',...
                ["energy_period"],[44.624076,-124.280097;43.489171,-125.152137],...
                2010,api_key);
            file = '../../examples/data/wave/hindcast/hindcast_1hr_data.csv';
            meta = '../../examples/data/wave/hindcast/hindcast_1hr_meta.csv';
            expected_data = readtable(file,'delimiter',',');
            expected_meta = readtable(meta);
            expected_data.time_index = datetime(expected_data.time_index,'InputFormat','yyyy-MM-dd HH:mm:ssXXX',...
                'TimeZone','UTC');

            assertEqual(testCase,expected_data.time_index,hindcast_data(1).time);
            assertEqual(testCase,expected_data.energy_period_0,hindcast_data(1).energy_period,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.latitude(1),hindcast_data(1).metadata.latitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.longitude(1),hindcast_data(1).metadata.longitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.water_depth(1),hindcast_data(1).metadata.water_depth,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.timezone(1),hindcast_data(1).metadata.timezone);
            assertEqual(testCase,expected_meta.jurisdiction{1},hindcast_data(1).metadata.jurisdiction);
            assertEqual(testCase,expected_meta.distance_to_shore(1),hindcast_data(1).metadata.distance_to_shore,'RelTol',0.000001);

            assertEqual(testCase,expected_data.time_index,hindcast_data(2).time);
            assertEqual(testCase,expected_data.energy_period_1,hindcast_data(2).energy_period,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.latitude(2),hindcast_data(2).metadata.latitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.longitude(2),hindcast_data(2).metadata.longitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.water_depth(2),hindcast_data(2).metadata.water_depth,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.timezone(2),hindcast_data(2).metadata.timezone);
            assertEqual(testCase,expected_meta.jurisdiction{2},hindcast_data(2).metadata.jurisdiction);
            assertEqual(testCase,expected_meta.distance_to_shore(2),hindcast_data(2).metadata.distance_to_shore,'RelTol',0.000001);
        end

        function test_WPTO_point_multiparm(testCase)

            assumeFail(testCase, "API key usage saturated - temporarily disabling")

            api_key = '3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf';
            hindcast_data = request_wpto('3-hour',...
                ["mean_absolute_period","significant_wave_height"],[44.624076,-124.280097],...
                1996,api_key);
            file = '../../examples/data/wave/hindcast/hindcast_3hr_data.csv';
            meta = '../../examples/data/wave/hindcast/hindcast_3hr_meta.csv';
            expected_data = readtable(file,'delimiter',',');
            expected_meta = readtable(meta);
            expected_data.time_index = datetime(expected_data.time_index,'InputFormat','yyyy-MM-dd HH:mm:ssXXX',...
                'TimeZone','UTC');

            assertEqual(testCase,expected_data.time_index,hindcast_data.time);
            assertEqual(testCase,expected_data.mean_absolute_period_0,hindcast_data.mean_absolute_period,'RelTol',0.000001);
            assertEqual(testCase,expected_data.significant_wave_height_0,hindcast_data.significant_wave_height,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.latitude(1),hindcast_data.metadata.latitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.longitude(1),hindcast_data.metadata.longitude,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.water_depth(1),hindcast_data.metadata.water_depth,'RelTol',0.000001);
            assertEqual(testCase,expected_meta.timezone(1),hindcast_data.metadata.timezone);
            assertEqual(testCase,expected_meta.jurisdiction{1},hindcast_data.metadata.jurisdiction);
            assertEqual(testCase,expected_meta.distance_to_shore(1),hindcast_data.metadata.distance_to_shore,'RelTol',0.000001);
        end

    end

end
