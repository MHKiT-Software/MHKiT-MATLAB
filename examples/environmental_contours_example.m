%% MHKiT Environmental Contours
% Environmental contours of extreme sea states can be used as a part of reliability-based
% design for offshore structures, including wave energy converters (WECs). Environmental
% contours provide estimations of extreme sea states based on short-term data
% (e.g., 10 years used to estimate a 100-year event). These environmental contours
% describe extreme sea states by characterizing the resource, defining sea states
% for extreme condition analysis, and developing a framework for analyzing survivability
% of a design.
%
% MHKiT includes functions adapted from the <https://github.com/WEC-Sim/WDRT
% WDRT> for creating environmental contours of extreme sea states using a principal
% component analysis (PCA) methodology, with additional improvements for characterizing
% the joint probability distribution of sea states. As a demonstration, this notebook
% will walk through the following steps to find a 100-year sea state for NDBC
% buoy 46022 using 16 years of spectral wave density data.
%%
% # Request Spectral Wave Density Data from NDBC
% # Calculate Hm0 and Te using the requested data
% # Find the data's 100-year contour
% # Plot the data and the 100-year contour

%% 1. Request Spectral Wave Density Data from NDBC
% MHKiT can be used to request historical data from the National Data Buoy Center
% (<https://www.ndbc.noaa.gov/ NDBC>). This process is split into the following
% steps:
%%
% * Query available NDBC data
% * Select years of interest
% * Request Data from NDBC
%
% *Query available NDBC data*
%
% The NDBC_available_data function requires a parameter to be specified and
% optionally the user may provide a station ID as a string. We are interested
% in historical spectral wave density data |'swden'| (from which we may calculate
% Hm0 and Te). Additionally, we will specify the buoy number as |'46022'| to only
% return data associated with this site.

% Specify the parameter as spectral wave density and the buoy number to be 46022
parameter = 'swden';
buoy_number = '46022';
available_data= NDBC_available_data(parameter,"buoy_number", buoy_number);
available_data
%%
% *Select years of interest*
%
% The |NDBC_available_data| function has returned a Table with columns 'Station_id',
% 'year', and 'file'. In this case, the years returned from NDBC|_available_data|
% span 1996 to the last complete year the buoy was operational (currently 2019
% for 46022). For demonstration, we have decided we are interested in the data
% between the years 1996 and 2012 so we will create a new |filenames_of_interest|
% variable which only contains filenames of years less than 2013.

% Slice the available data to only include through year 2012
rows = (available_data.year < 2013) ;
filenames_of_interest = available_data.file(rows);
filenames_of_interest

%%
% *Request Data from NDBC*
%
% To get the NDBC data we can use the NDBC|_request_data| function to iterate
% over each buoy id and year in passed filenames. This function will return the
% parameter data as a structure of structures which may be accessed by buoy id
% and then the year for multiple buoys or just the year for a single buoy. An
% additional data column called 'time' is created with time in datetime format.

ndbc_requested_data = NDBC_request_data(parameter, filenames_of_interest);

%% 2. Calculate Hm0 and Te using the NDBC Data
%
% A sea state may be characterized by significant wave height (Hm0) and energy
% period (Te). Using the historical spectral wave density data from NDBC, we can
% calculate these variables using MHKiT. Both Hm0 and Te return a single value
% for a given time (e.g., DateTime index).

Hm0 = [];
Te = [];
for field = fieldnames(ndbc_requested_data)'
    Hm0 = [Hm0 ; significant_wave_height(ndbc_requested_data.(field{1}))];
    Te = [Te ; energy_period(ndbc_requested_data.(field{1}))];
end

%% 3. Find the 100-year contour line
%
% With the sea state data calculated, we can now use the modified I-FORM method
% to define reliability for a 100-year sea state based on the 17 years of spectral
% wave density data obtained from NDBC for buoy 46022. Reliability is the likelihood
% that a certain event will not occur in a given period. The period will define
% a line of constant probability in the joint probability of Hm0 and Te but individually
% each component different reliability (marginal distribution) which we can find
% by evaluating a normal cumulative distribution function (CDF). This CDF returns
% each component's quantiles along the iso-reliability line that finally allows
% us to calculate each sea state value (e.g., the 100-year contour values for
% Hm0 and Te).
%
% For more detail on the environmental contour method used here please refer
% to: <https://www.sciencedirect.com/science/article/abs/pii/S0029801815006721
% Eckert-Gallup et. al 2016>
%
% To apply the environmental contours function we will specify a 100-year sea
% state, the sea state data (Hm0, Te), and the time difference between measurements
% (dt in seconds).

% Return period (years) of interest
period = 100 ;

% Remove Hm0 Outliers and NaNs
filter = Hm0 < 20;
Hm0 = Hm0(filter);
Te = Te(filter);
[row, col] = find(~isnan(Te));
Hm0 = Hm0(row);
Te = Te(row);
[row, col] = find(~isnan(Hm0));
Hm0 = Hm0(row);
Te = Te(row);

% Delta time of sea-states in seconds
dt = ndbc_requested_data.year_1996.time(2)- ndbc_requested_data.year_1996.time(1);
dt= seconds(dt);

% Get the contour values
contour = environmental_contours(Hm0, Te, dt, period, 'PCA');

%% 4. Plot overlay of the data and contour
%
% Lastly we can use the MHKiT graphics module to create a contour plot which
% shows the data and resultant contour line.

figure('Position', [100, 100, 1600, 600]);

plot_environmental_contours(Te,Hm0,contour.contour2,contour.contour1,"x_label",...
    'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
    "contour_label",'100 Year Contour');

%% Other Contour Methods
%
% MHKiT has parametric, nonparametric, and Kernel Density Estimation  methods
% for calculating environmental contours housed within the |resource.environmental_contours|
% function. We can compare other copulas to our PCA method simply by adding more
% methods. A single string method can be applied if only one  copula is of interest
% as was done above for the PCA method or multiple  methods can be sent in using
% a list of strings. If multiple methods of the same type are desired is recommended
% to run them at the same time if possible as it will reduce the computational
% expense by utilizing the  common computational calls across the copulas. In
% the example below we  will compare the parametric and nonparametric Gaussain
% contours to the PCA method we ran previously.
%
% We start by using the default settings for our Gaussian methods. The environmental_contours
% function returns a dictionary with component 'x1' and 'x2' for each method.
% E.g. if the environmental_contours function were called with the 'gaussian'
% method then environmental_contours function would return a dictionary with keys
% ['gaussian_x1', 'gaussian_x2']. The copula methods are a generalized mathematical
% method and therefore 'x1' and 'x2' are used in place of Hm0 and Te for the component
% values. 'x1' refers to the first array passed and 'x2' refers to the second
% array passed. In the example below 'x1' would refer to the Hm0 component of
% the coupla and 'x2' would refer to Te.

gauss = environmental_contours(Hm0, Te, dt, period, "gaussian");
nongauss = environmental_contours(Hm0, Te, dt, period, "nonparametric_gaussian");

Tes = [contour.contour2; gauss.contour2; nongauss.contour2];
Hm0s = [contour.contour1; gauss.contour1; nongauss.contour1];

plot_environmental_contours(Te,Hm0,Tes,Hm0s,"x_label",...
    'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
    "contour_label",["PCA","Gaussian","Nonparametric Gaussian"]);
