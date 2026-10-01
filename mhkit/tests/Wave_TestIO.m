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

            assumeFail(testCase, "Hindcast tests skipped - slow network calls")
            api_key = '3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf';
            hindcast_data = request_wpto('3-hour',...
                ["energy_period"],[44.624076,-124.280097;43.489171,-125.152137],...
                2010,api_key);
            file = '../../examples/data/wave/hindcast/hindcast_3hr_multiloc_data.csv';
            meta = '../../examples/data/wave/hindcast/hindcast_3hr_multiloc_meta.csv';
            expected_data = readtable(file,'delimiter',',');
            expected_meta = readtable(meta);
            expected_data.time_index = datetime(expected_data.time_index,'InputFormat','yyyy-MM-dd HH:mm:ssXXX',...
                'TimeZone','UTC');

            assertEqual(testCase,expected_data.time_index,hindcast_data(1).time);
            assertEqual(testCase,expected_data.energy_period_0,hindcast_data(1).energy_period,'RelTol',0.000001);
            verify_site_metadata(testCase, expected_meta, 1, hindcast_data(1).metadata);

            assertEqual(testCase,expected_data.time_index,hindcast_data(2).time);
            assertEqual(testCase,expected_data.energy_period_1,hindcast_data(2).energy_period,'RelTol',0.000001);
            verify_site_metadata(testCase, expected_meta, 2, hindcast_data(2).metadata);
        end

        function test_WPTO_omnidirectional(testCase)

            assumeFail(testCase, "Hindcast tests skipped - slow network calls")

            api_key = '3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf';
            hindcast_data = request_wpto('3-hour',...
                ["omni-directional_wave_power"],[44.624076,-124.280097;43.489171,-125.152137],...
                2010,api_key);

            % Verify hindcast data is not empty
            testCase.verifyNotEmpty(hindcast_data, 'Hindcast data should not be empty.');

            expected_point_1_file = "../../examples/data/wave/hindcast/hindcast_test_omni_point_1.csv";
            expected_point_1 = readtable(expected_point_1_file,'delimiter',',');

            expected_point_2_file = "../../examples/data/wave/hindcast/hindcast_test_omni_point_2.csv";
            expected_point_2 = readtable(expected_point_2_file,'delimiter',',');

            % Convert missing csv timezone to UTC
            expected_point_1.time = datetime(expected_point_1.time, 'TimeZone', 'UTC');
            expected_point_2.time = datetime(expected_point_2.time, 'TimeZone', 'UTC');

            % Verify downloaded omni-directional_wave_power matches the expected values
            assertEqual(testCase, expected_point_1.omni_directional_wave_power,hindcast_data(1).omni_directional_wave_power, 'RelTol',0.000001);
            assertEqual(testCase, expected_point_2.omni_directional_wave_power,hindcast_data(2).omni_directional_wave_power, 'RelTol',0.000001);

            % Verify time matches the expected values
            assertEqual(testCase, expected_point_1.time,hindcast_data(1).time);
            assertEqual(testCase, expected_point_2.time,hindcast_data(2).time);
        end

        function test_WPTO_point_multiparm(testCase)

            assumeFail(testCase, "Hindcast tests skipped - slow network calls")

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
            verify_site_metadata(testCase, expected_meta, 1, hindcast_data.metadata);
        end

    end

end

function verify_site_metadata(testCase, expected_meta, row, metadata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare request_wpto site metadata with the expected values
%
% The expected metadata was read from the hindcast source files on AWS
% S3. As of 2026-09-30, HSDS returns an empty response for the
% West_Coast 3-hour meta dataset, so request_wpto reads latitude,
% longitude, and water_depth from other datasets and leaves
% distance_to_shore, timezone, and jurisdiction missing. Both cases are
% checked so the tests stay valid when HSDS is fixed.
%
% Parameters
% ------------
% testCase : matlab.unittest.TestCase
%   Test case
% expected_meta : table
%   Expected metadata, one row per site
% row : double
%   Row of expected_meta to compare
% metadata : struct
%   Site metadata returned by request_wpto
%
% Returns
% ---------
% None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assertEqual(testCase, expected_meta.latitude(row), metadata.latitude, 'RelTol', 0.000001);
    assertEqual(testCase, expected_meta.longitude(row), metadata.longitude, 'RelTol', 0.000001);
    assertEqual(testCase, expected_meta.water_depth(row), metadata.water_depth, 'RelTol', 0.000001);
    if isnan(metadata.timezone)
        % meta dataset unavailable from HSDS
        assertTrue(testCase, isnan(metadata.distance_to_shore));
        assertEqual(testCase, metadata.jurisdiction, '');
    else
        assertEqual(testCase, expected_meta.timezone(row), metadata.timezone);
        assertEqual(testCase, expected_meta.jurisdiction{row}, metadata.jurisdiction);
        assertEqual(testCase, expected_meta.distance_to_shore(row), metadata.distance_to_shore, 'RelTol', 0.000001);
    end
end
