%% Extreme Conditions Modeling - Contour Approach
%
% Extreme conditions modeling consists of identifying the expected  extreme
% (e.g. 100-year) response of some quantity of interest, such as  WEC motions
% or mooring loads. Three different methods of estimating extreme conditions were
% adapted  from <https://github.com/WEC-Sim/WDRT WDRT>: full sea state approach,
% contour approach, and MLER design wave. This noteboook presents the contour
% approach.
%
% The contour approach consists of the following steps:
%%
% # Sample the environmental contour for the return period of interest at N
% points near the maximum _*H*_.
% # For each sample _*(H, T)*_ calculate the short-term (e.g. 3-hours) extreme
% for the quantity of interest (e.g. WEC motions or mooring tension).
% # Select the sea state with the highest short-term expected value for the
% quantity of interest as the design sea state, and some metric (e.g. 95th percentile)
% of the distribution as the extreme/design response.
%%
% *NOTE:* Prior to running this example it is recommended to become familiar
% with |environmental_contours_example.ipynb| and |short_term_extremes_example.ipynb|
% since some code blocks are adapted from those examples and used here without
% the additional description.
%
%% Obtain and Process NDBC Buoy Data
% The first step will be obtaining the environmental data and creating the contours.
% See |environmental_contours_example.ipynb| for more details and explanations
% of how this is being done in the following code block.

% Specify the parameter as spectral wave density and the buoy number to be 46022
parameter = 'swden';
buoy_number = '46022';
available_data= NDBC_available_data(parameter,"buoy_number", buoy_number);

% Slice the available data to only include through year 2012
rows = (available_data.year < 2013) ;
filenames_of_interest = available_data.file(rows);

ndbc_requested_data = NDBC_request_data(parameter, filenames_of_interest);

Hm0 = [];
Te = [];
for field = fieldnames(ndbc_requested_data)'
    Hm0 = [Hm0 ; significant_wave_height(ndbc_requested_data.(field{1}))];
    Te = [Te ; energy_period(ndbc_requested_data.(field{1}))];
end

% Remove Hm0 Outliers and NaNs
filter = Hm0 < 20;
Hm0 = Hm0(filter);
Te = Te(filter);
[row, ~] = find(~isnan(Te));
Hm0 = Hm0(row);
Te = Te(row);
[row, col] = find(~isnan(Hm0));
Hm0 = Hm0(row);
Te = Te(row);

% Delta time of sea-states in seconds
dt = ndbc_requested_data.year_1996.time(2)- ndbc_requested_data.year_1996.time(1);
dt = seconds(dt);

%% 1. Sampling
%
% The first step is to create the environmental contour for the return period
% of interest, which is 100-years in this case, using the |wave.contours.environmental_contours|.
% See  |environmental_contours_example.ipynb| for more details on using this function.
%
% Next we choose 5 sample sea states along this contour. We choose 5  equally
% spaced energy period values, between 15-22 seconds, to sample  the contour.
% We then obtain the significant wave height for each sample  using the |wave.contours.samples_contour|
% function.

% Get the contour values
period = 100; % 100-year contour
contour = environmental_contours(Hm0, Te, dt, period, "PCA");
hs_contour = contour.contour1;
te_contour = contour.contour2;
% 5 samples
te_samples = linspace(15, 22, 5);
hs_samples = samples_contour(te_samples, te_contour, hs_contour);
%%
% We now plot the contour and samples.

plot_environmental_contours(Te,Hm0,te_contour,hs_contour,"x_label",...
    'Energy Period (s)', "y_label",'Significant Wave Height (m)',"data_label",'NDBC 46022',...
    "contour_label",'100 Year Contour');
hold on;
scatter(te_samples, hs_samples,"red","filled","o"); legend("data","100-year contour","samples")
hold off;

%% 2. Short-Term Extreme Distributions
%
% Many different methods for short-term extremes were adapted from WDRT, and
% a summary and examples can be found in |short_term_extremes_example.ipynb|.
% The response quantity of interest is typically related to the WEC  itself, e.g.
% maximum heave displacement, PTO extension, or load on the  mooring lines. This
% requires running a simulation (e.g. WEC-Sim) for each of the 5 sampled sea states
% _*(H, T)*_. For the sake of example we will consider the wave elevation as the
% quantity of interest (can be thought as a proxy for heave motion in this example).
% Wave elevation time-series for a specific sea state can be created  quickly
% without running any external software.
%
% *NOTE:* The majority of the for-loop below is simply  creating the synthetic
% data (wave elevation time series). In a realistic case the variables |time|
% and |data| describing  each time series would be obtained externally, e.g. through
% simulation  software such as WEC-Sim or CFD. For this reason the details of
% creating the synthetic data are not presented here, instead assume for each
% sea  state there is time-series data available.
%
% The last lines of the for-loop create the short-term extreme distribution
% from the time-series using the |loads.extreme.short_term_extreme| function.
% The short-term period will be 3-hours and we will use 1-hour  "simulations"
% and the Weibul-tail-fitting method to estimate the 3-hour  short-term extreme
% distributions for each of the 5 samples.
%
% For more details see |short_term_extremes_example.ipynb| and
%
% _[3] Michelén Ströfer, Carlos A., and Ryan Coe. 2015. “Comparison of  Methods
% for Estimating Short-Term Extreme Response of Wave Energy  Converters.” In OCEANS
% 2015 - MTS/IEEE Washington, 1–6. IEEE._

% create the short-term extreme distribution for each sample sea state
t_st = 3.0 * 60.0 * 60.0;
gamma = 3.3;
t_sim = 1.0 * 60.0 * 60.0;

n = length(hs_samples);

for i=1:length(hs_samples)
    tp = te_samples(i) / (0.8255 + 0.03852 * gamma - 0.005537 * gamma^2 + 0.0003154 * gamma^3);
    % time & frequency arrays
    df = 1.0 / t_sim;
    T_min = tp / 10.0;  % s
    f_max = 1.0 / T_min;
    Nf = int32(f_max / df) + 1;
    time = linspace(0, t_sim, 2 * Nf + 1);
    f = linspace(0.0, f_max, Nf);
    fprintf("Sea state %g/%g. (Hs, Te) = (%.2f m, %.2f s). Tp = %.2f s\n",i,n,hs_samples(i),te_samples(i),tp);
    % spectrum
    S = jonswap_spectrum(f, tp, hs_samples(i), gamma);
    % 1-hour elevation time-series
    data = surface_elevation(S, time);
    % 3-hour extreme distribution
    ste_all(i) = short_term_extreme(time, data.elevation, t_st, "peaks_weibull_tail_fit", 1, "expect");
    ste_95(i) = short_term_extreme(time, data.elevation, t_st, "peaks_weibull_tail_fit", 0.95, "ppf");
end

%% 3. Select Design Sea State and Design Response
%
% Finally, we will choose the sea state with the largest expected value as the
% design sea state. In the code snippet above, we stored the expected value for
% each sea state in |ste_all|. We will then use the 95th percentile of the short-term
% extreme distribution as the extreme design response.
%
% First we find the sampled sea state with the largest expected value.

% find max expected wave height
[~,max_ind] = max(ste_all);
hs_design = hs_samples(max_ind);
te_design = te_samples(max_ind);
fprintf("Design sea state (Hs, Te): (%.2f m, %.2f s)", hs_design, te_design)

%%
% Next, we choose the 95th percentile of the short-term extreme distribution
% at this sea state as the extreme design response.

response_design = ste_95(max_ind);
fprintf("Design (extreme) response: %.2f meters", response_design)