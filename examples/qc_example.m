%% Example: MHKiT-MATLAB Quality Control Module
%
% The following example runs a simple quality control analysis on wave elevation
% data using the <https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.qc.html
% MHKiT QC module>. The data file used in this example is stored in the <https://github.com/MHKiT-Software/MHKiT-MATLAB/tree/master/examples/data/qc
% /mhkit/examples/data/qc> directory.
%
%% Import Data
%
% Load the data into a table. The data includes several issues, including timestamps
% that are out of order, corrupt data with values of -999, data outside expected
% range, and stagnant data.

data = readtable('./examples/data/qc/wave_elevation_data.csv');
disp(data)
set(gcf, 'Position',  [100, 100, 2000, 600]); % setting the size of the plots
plot(data.Time,[data.probe1,data.probe2,data.probe3]);
ylim([-60 60]);
xlabel('Time')
legend('probe 1', 'probe 2', 'probe 3');

%%
% The MHKiT-MATLAB QC module utilizes structures. We need to convert the table
% of data into a structure.

dataset.time=data.Time;
dataset.values = table2array(data(:,2:4));
dataset

%% Quality Control Tests
%
% The following quality control tests are used to identify timestamp issues,
% corrupt data, data outside expected range, and stagnant data.
%
% Each quality control tests results in the following information in a results
% structure:
%%
% * Cleaned data that has _NaN_ in place of data that did not pass the quality control test
% * Boolean mask with True/False that indicates if each data point passed the quality control test

%%
% *Check Timestamps*
%
% Quality control analysis generally starts by checking the timestamp index
% of the data.
%
% The following test checks to see if:
%
% # the data contains duplicate timestamps,
% # timestamps are not monotonically increasing
% # timestamps occur at irregular intervals (an interval of 0.002s is expected for this data).
%
% If duplicate timestamps are found, the resulting structure keeps the first
% occurrence. If timestamps are not monotonic, the timestamps in the resulting
% structure are reordered.

% define expected frequency of the data, in seconds
frequency = 0.002;

% run the timestamp quality control test
results = check_timestamp(dataset,frequency);

%%
% The cleaned data and boolean mask are shown below.

disp(results);
plot(results.time,results.values);
ylim([-60 60]);
xlabel('Time')
legend('probe 1', 'probe 2', 'probe 3');

%%
% *Check for Corrupt Data*
%
% In the following quality control tests, the cleaned data from the previous
% test are used as input to the subsequent test. For each quality control test,
% a plot of the cleaned data is shown.
%
% The quality control test below checks for corrupt data, indicated by a value
% of -999.

%define corrupt values
corrupt_values = {-999};

% run the corrupt data QC test
results = check_corrupt(results,corrupt_values);

% Plot the cleaned data
plot(results.time,results.values);
ylim([-60 60]);
xlabel('Time')
legend('probe 1', 'probe 2', 'probe 3');

%%
% *Check for Data Outside Expected Range*
%
% The next quality control test checks for data that is greater than 50 or less
% than -50.  Note that expected range tests can also be used to compare measured
% values to a model, or analyze the expected relationships between data columns.

% define expected lower and upper bound
expected_bounds = [-50,50];

% run expected range QC test
results = check_range(results,expected_bounds);

% Plot the cleaned data
plot(results.time,results.values);
ylim([-60 60]);
xlabel('Time')
legend('probe 1', 'probe 2', 'probe 3');

%%
% *Check for stagnant data*
%
% The final quality control test checks for stagnant data by looking for data
% that changes by less than 0.001 within a 0.02 second moving window.

% define the lower bound (no upper bound is specified in this example)
expected_bound = {0.001,py.None};

% Define the moving window, in seconds
window = 0.02;

% input the data, boundaries, py.None indicates apply to all columns, and window
results = check_delta(results,expected_bound,window);

% plot cleaned data
plot(results.time,results.values);
ylim([-60 60]);
xlabel('Time')
legend('probe 1', 'probe 2', 'probe 3');

%% Cleaned Data
%
% The cleaned data can be used directly in MHKiT analysis, or the  missing values
% can be replaced using various methods before analysis is  run. Data replacement
% strategies are generally defined on a case-by-case basis. MATLAB includes methods
% to interpolate, replace, and fill missing values.

% Extract final cleaned data for MHKiT analysis
results
disp(results.values)
