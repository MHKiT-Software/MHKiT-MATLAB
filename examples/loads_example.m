%% Example: MHKiT-MATLAB Loads Example
%
% The following example will help familiarize you with some of the functions
% in the <https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.loads.html MHKiT
% Loads Module> that you can use to assist you in your loads analysis.

%% Import Loads Data
%
% In this example, since there is no known field loads data for MHK devices,
% this data comes from a land-based wind turbine. We have a database of 331 files,
% each containing 10 minutes of data sampled at 50Hz.
%
% As a start, let's look at the data for one of these files to figure out what
% formatting we need to apply.

data = readtable('./examples/data/loads/data_loads_example.csv');
disp(data)

%%
% *Format Loads Data with datetime index*
%
% Our file above shows us that we have two references to time, but neither are
% in the right format. The "Timestamp" column is what will give us the datetime
% index that we are looking for, but we first need to convert it from Microsoft
% Excel format to unix time.

newtime = excel_to_datetime(data.Timestamp);
data.time = newtime

%% Loads Analysis
%
% Now that we have our loads data in the correct format for MHKiT, we do some
% analysis.

%%
% *Damage Equivalent Loads*
%
% Let's say that we wanted to investigate fatigue. We can do this by calculating
% short-term damage equivalent loads (DELs). To use this function, we need the
% variable we want to analyze and its corresponding material slope factor. In
% this instance, we want to just look at our tower base moment and our blade 1
% root flap moment. Our tower is steel while our blade is composite so they will
% have different material slopes.
%
% We call our function and apply the default inputs of using at least 100 bins
% for the load ranges and we let t=600 seconds so that we get an equivalent 1Hz
% DEL for our 10-minute file.

% use the fatigue function to calculate the damage equivalent load for this file
DEL_tower = damage_equivalent_load(data.TB_ForeAft,4,'bin_num',100,'data_length',600)
DEL_blade = damage_equivalent_load(data.BL1_FlapMom,10,'bin_num',100,'data_length',600)

%%
% *Calculate Loads Statistics*
%
% Another important part of loads analysis is looking at statistics. Here, we
% use another function to help us calculate the mean, max, min, and std for this
% 10-minute file. Per standards, a valid statistical window has to be consecutive
% in time with the correct number of datapoints. So, if this 10 minute file did
% not meet this criteria, then no stats would be generated and a warning message
% would appear.
%
% NOTE: Sometimes individual files may contain enough data for multiple statistical
% windows. This function can still handle this scenario as long as the correct
% inputs are specified.

% First we need to convert the table to a structure
data = removevars(data,{'Timestamp','Time'});
datast = table2struct(data,'ToScalar',true)
% calculate the means, maxs, mins, and stdevs of this file
stats = get_statistics(datast,50);
% show the mean results
disp(stats.mean)

%%
% At this point, it would be nice to start visualizing some of this data. In
% order to do this, we need to calculate the stats and DELs for all the files
% in our database. In this case, it would be done through a loop that imports
% each file and applies all the functions we just saw.  To speed things up, this
% was already done so we just need to import the results.

% load statistics results from csv
dmeans = table2struct(readtable('./examples/data/loads/data_loads_means.csv'),'ToScalar',true)
dmaxs = table2struct(readtable('./examples/data/loads/data_loads_maxs.csv'),'ToScalar',true);
dmins = table2struct(readtable('./examples/data/loads/data_loads_mins.csv'),'ToScalar',true);
dstd = table2struct(readtable('./examples/data/loads/data_loads_std.csv'),'ToScalar',true);

%%
% *Visualizing Load Statistics*
%
% Now that we have all of our stats, let's make a scatter plot that can give
% us a visual. Using the |plot_statistics| function, we can quickly create a standard
% scatter plot showing how load variables trend with wind speed.  Using this we
% can quickly spot expected trends and track down outliers.

figure = plot_statistics(dmeans.uWind_80m,dmeans.BL1_FlapMom,dmins.BL1_FlapMom,dmaxs.BL1_FlapMom,"y_stdev",dstd.BL1_FlapMom,"xlabel",'Wind Speed [m/s]',"ylabel",'Blade Flap Mom [kNm]');
figure2 = plot_statistics(dmeans.uWind_80m,dmeans.TB_ForeAft,dmins.TB_ForeAft,dmaxs.TB_ForeAft,"y_stdev",dstd.TB_ForeAft,"xlabel",'Wind Speed [m/s]',"ylabel",'Tower Base Mom [kNm]');

%%
% Another common step is to bin the statistical data. This can easily be done
% with the |bin_statistics| function from the loads module shown below. A warning message
% will show if there are any bins that were not filled.

% create array containing wind speeds to use as bin edges
b_edges = 3:1:25 ;
% apply function for means, maxs, mins, and DELs
wind_means = bin_statistics(dmeans,dmeans.uWind_80m,b_edges);
wind_max = bin_statistics(dmaxs,dmeans.uWind_80m,b_edges);
wind_min = bin_statistics(dmins,dmeans.uWind_80m,b_edges);
wind_DEL = bin_statistics(dmaxs,dmeans.uWind_80m,b_edges);

wind_min.averages

%%
% Now let's make some more plots with the binned data. Here we use the binned
% data and corresponding standard deviations as inputs to the |plot_bin_statistics| function.

bcenters = 3.5:1:24.5
plot_bin_statistics(bcenters,wind_means.averages.TB_ForeAft,wind_max.averages.TB_ForeAft,wind_min.averages.TB_ForeAft,wind_means.std.TB_ForeAft,wind_max.std.TB_ForeAft,wind_min.std.TB_ForeAft,"xlabel",'Wind Speed [m/s]',"ylabel",'TB-ForeAft',"title",'Binned Stats');
