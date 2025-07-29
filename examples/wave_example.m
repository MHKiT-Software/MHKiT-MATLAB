%% Example: MHKiT-MATLAB Wave Module
% The following example runs an application of the <https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.wave.html
% MHKiT wave module> to 1) generate a capture length matrix, 2) calculate MAEP,
% and 3) plot the scatter diagrams.
%% Load NDBC Data
% We can use MHKiT to load data downloaded from <https://www.ndbc.noaa.gov/
% https://www.ndbc.noaa.gov>.

relative_file_name = './examples/data/wave/data.txt';
current_dir = fileparts(matlab.desktop.editor.getActiveFilename);
full_file_name = fullfile(current_dir, relative_file_name);
ndbc_data = read_NDBC_file(full_file_name);
disp(ndbc_data)

%% Compute Wave Metrics
% We will now use MHKiT to compute the significant wave height, energy period,
% and energy flux.

% Compute the energy periods from the NDBC spectra data
Te = energy_period(ndbc_data);
disp(Te)
% Compute the significant wave height from the NDBC spectra data
Hm0 = significant_wave_height(ndbc_data);
disp(Hm0)
% Set water depth to 60 m
h = 60;

% Compute the energy flux  from the NDBC spectra data and water depth
J = energy_flux(ndbc_data,h);
disp(J)
%% Generate Random Power Data
% For demonstration purposes, this example uses synthetic power data generated
% from statistical distributions. In a real application, the user would provide
% power values from a WEC.

% generating 1,000,000 random power values
Power = randi([40,200], 743,1);
disp(Power)
%% Capture Length Matrices
% The following operations create capture length matrices, as specified by the
% IEC/TS 62600-100. But first, we need to calculate capture length and define
% bin centers. The mean capture length matrix is printed below. Keep in mind that
% this data has been artificially generated, so it may not be representative of
% what a real-world scatter diagram would look like.

% calculating capture length with power and wave flux in vectors
L = capture_length(Power,J)

% Need to set our Hm0 and Te bins for the capture length matrix
Hm0_bins = -0.5:0.5:max(fix(Hm0))+0.5; % Input is min, max, and n indices for vector
Hm0_bins = Hm0_bins+0.25 ;
Te_bins = 0:1:max(fix(Te));
Te_bins = Te_bins+0.5;

% Calculate the necessary capture length matrices for each statistic based
% on IEC/TS 62600-100
clmat.mean = capture_length_matrix(Hm0,Te,L,"mean",Hm0_bins,Te_bins);
clmat.std = capture_length_matrix(Hm0,Te,L,"std",Hm0_bins,Te_bins);
clmat.count = capture_length_matrix(Hm0,Te,L,"count",Hm0_bins,Te_bins);
clmat.min = capture_length_matrix(Hm0,Te,L,"min",Hm0_bins,Te_bins);
clmat.max = capture_length_matrix(Hm0,Te,L,"max",Hm0_bins,Te_bins);

% Calculate the frequency matrix for convenience
clmat.freq = capture_length_matrix(Hm0,Te,L,"frequency",Hm0_bins,Te_bins);
%%
% Let's see what the data in the mean matrix looks like. Keep in mind that this
% data has been artificially generated, so it may not be representative of what
% a real-world scatter diagram would look like.

disp(clmat.mean.values)
%% Power Matrices
% As specified in IEC/TS 62600-100, the power matrix is generated from the capture
% length matrix and wave energy flux matrix, as shown below

% Create wave energy flux matrix using mean
jmat = wave_energy_flux_matrix(Hm0,Te,J,"mean",Hm0_bins,Te_bins);

% Create power matrix using mean
avg_power_mat = power_matrix(clmat.mean, jmat);

% Create power matrix using standard deviation
std_power_mat = power_matrix(clmat.std, jmat);
%%
% The |capture_length_matrix| function can also be used as an arbitrary scatter
% plot generator. To do this, simply pass a different array in the place of capture
% length (L). For example, while not specified by the IEC standards, if the user
% doesn't have the omnidirectional wave flux, the average power matrix could hypothetically
% be generated in the following manner:

avgpowmat_not_standard = capture_length_matrix(Hm0,Te,Power,'mean',Hm0_bins,Te_bins);
%% MAEP
% There are two ways to calculate mean annual energy production (MAEP). One
% is from capture length and wave energy flux matrices, the other is from time
% series data, as shown below.
%
%

% Calculate maep from timeseries
maep_timeseries = mean_annual_energy_production_timeseries(L,J)
% Calculate maep from matrix
maep_matrix = mean_annual_energy_production_matrix(clmat.mean, jmat , clmat.freq)
%% Graphics
% The graphics function |plot_matrix| can be used to visualize results. It is
% important to note that the plotting function assumes the step size between bins
% to be linear.

% Plot the capture length matrix
p1 = plot_matrix(clmat.mean,"Capture Length", "annotate", false);
