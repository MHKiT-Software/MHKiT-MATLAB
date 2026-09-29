%% Example: MHKiT-MATLAB Wave Module
% The following example runs an application of the <https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.wave.html
% MHKiT wave module> to 1) generate a capture width matrix, 2) calculate MAEP,
% and 3) plot the scatter diagrams.
%% Load CDIP Wave Measurement Data PacWave North, CDIP 277
% This example uses one year (2025) of wave spectra from the PacWave
% North buoy off Newport, Oregon, published by CDIP as
% <https://cdip.ucsd.edu/m/products/?stn=277p1 station 277>. The spectra are
% requested from the CDIP THREDDS server with |cdip_request_parse_workflow|.
%
% A saved copy of that request is shipped with the examples so this script
% runs without network access. A working request is shown commented out
% below.

station_number = '277';

% cdip = cdip_request_parse_workflow('station_number', station_number, ...
%     'data_type', 'historic', ...
%     'start_date', '2025-01-01', 'end_date', '2025-12-31', ...
%     'parameters', {'waveEnergyDensity', 'metaWaterDepth'});
% % waveEnergyDensity is stored [time x frequency]; MHKiT uses frequency
% % down the rows and one spectrum per column
% S.frequency = double(cdip.metadata.wave.waveFrequency(:));
% S.spectrum = double(cdip.data.wave2D.waveEnergyDensity');
% S.time = cdip.data.wave.waveTime(:);
% S.time.TimeZone = 'UTC';
% % Water depth at the buoy [m], a scalar metadata variable in the NetCDF
% h = double(cdip.metadata.meta.metaWaterDepth);

% Saved copy of the request above
saved = load('./data/wave/cdip_277_2025.mat');
S.frequency = saved.cdip_277.frequency;
S.spectrum = double(saved.cdip_277.spectrum);  % stored as single, as CDIP publishes it
S.time = saved.cdip_277.time;
h = saved.cdip_277.water_depth;  % [m], from metaWaterDepth
fprintf('Using saved CDIP 277 data downloaded %s UTC\n', saved.cdip_277.downloaded);
disp(S)
time = S.time;

% Shared plot title suffix
title_suffix = sprintf('CDIP %s, %s', station_number, datetime(time(1), 'Format', 'yyyy'));
%% Compute Wave Metrics
% We will now use MHKiT to compute the significant wave height, energy period,
% and energy flux from each wave spectrum.

Hm0 = significant_wave_height(S);  % [m]
Te = energy_period(S);             % [s]
J = energy_flux(S, h);             % [W/m], uses the water depth

figure('Position', [100 100 900 700]);
tiledlayout(3, 1, 'TileSpacing', 'compact');
sgtitle(sprintf('Wave Metrics | %s', title_suffix))
nexttile; plot(time, Hm0); ylabel('H_{m0} [m]'); title('Significant Wave Height'); grid on
nexttile; plot(time, Te); ylabel('T_e [sec]'); title('Energy Period'); grid on
nexttile; plot(time, J / 1000); ylabel('J [kW/m]'); title('Wave Energy Flux'); grid on; xlabel('Time [UTC]')
%% Generate Random Power Data
% For demonstration purposes, this example uses synthetic power data generated
% from a uniform distribution. In a real application, the user would provide
% power values measured from a WEC.

rng(1);  % fixed seed so the example is repeatable
Power = randi([40, 200], numel(time), 1);  % [W]
%% Capture Width Matrices
% The following operations create capture width matrices, as specified by the
% IEC/TS 62600-100. But first, we need to calculate capture width and define
% bin centers. Keep in mind that the power has been artificially generated, so
% the capture width is not representative of a real WEC.

% calculating capture width with power and wave flux in vectors
CW = capture_width(Power, J);

% Bin centers with IEC TS 62600-100 widths of 0.5 m for Hm0 and 1 s for Te,
% so the edges fall on 0, 0.5, 1, ... m and 0, 1, 2, ... s.
Hm0_bins = 0.25:0.5:ceil(max(Hm0) / 0.5) * 0.5;
Te_bins = 0.5:1:ceil(max(Te));

% Calculate the necessary capture width matrices for each statistic based
% on IEC/TS 62600-100
cwmat.mean = capture_width_matrix(Hm0, Te, CW, "mean", Hm0_bins, Te_bins);
cwmat.std = capture_width_matrix(Hm0, Te, CW, "std", Hm0_bins, Te_bins);
cwmat.count = capture_width_matrix(Hm0, Te, CW, "count", Hm0_bins, Te_bins);
cwmat.min = capture_width_matrix(Hm0, Te, CW, "min", Hm0_bins, Te_bins);
cwmat.max = capture_width_matrix(Hm0, Te, CW, "max", Hm0_bins, Te_bins);

% Calculate the frequency matrix for convenience
cwmat.freq = capture_width_matrix(Hm0, Te, CW, "frequency", Hm0_bins, Te_bins);
%%
% Let's see what the data in the mean matrix looks like. Te bins run across
% the columns and Hm0 bins down the rows.

disp(cwmat.mean.values)
%% Power Matrices
% As specified in IEC/TS 62600-100, the power matrix is generated from the capture
% width matrix and wave energy flux matrix, as shown below

% Create wave energy flux matrix using mean
jmat = wave_energy_flux_matrix(Hm0, Te, J, "mean", Hm0_bins, Te_bins);

% Create power matrix using mean
avg_power_mat = power_matrix(cwmat.mean, jmat);

% Create power matrix using standard deviation
std_power_mat = power_matrix(cwmat.std, jmat);
%%
% The |capture_width_matrix| function can also be used as an arbitrary scatter
% plot generator. To do this, simply pass a different array in the place of capture
% width (CW). For example, while not specified by the IEC standards, if the user
% doesn't have the omnidirectional wave flux, the average power matrix could hypothetically
% be generated in the following manner:

avgpowmat_not_standard = capture_width_matrix(Hm0, Te, Power, 'mean', Hm0_bins, Te_bins);
%% MAEP
% There are two ways to calculate mean annual energy production (MAEP). One
% is from capture width and wave energy flux matrices, the other is from time
% series data, as shown below.

% Calculate maep from timeseries
maep_timeseries = mean_annual_energy_production_timeseries(CW, J)
% Calculate maep from matrix
maep_matrix = mean_annual_energy_production_matrix(cwmat.mean, jmat, cwmat.freq)
%% Graphics
% The graphics function |plot_matrix| can be used to visualize results. Each
% bin is drawn between its edges with the tick marks at the bin edges, and
% |trim_to_data| limits the axes to the bins that hold data.

% Plot the capture width matrix
figure('Position', [100 100 1000 800]);
ax1 = plot_matrix(cwmat.mean, "Capture Width", "zlabel", "Capture Width [m]", ...
    "trim_to_data", true, "font_size", 8);
title(ax1, sprintf('Mean Capture Width Matrix | %s', title_suffix));
%%
% The frequency matrix is the joint probability distribution (JPD) of the sea
% states: the fraction of records that fall in each Hm0 and Te bin. It is
% plotted here as a percentage with the empty bins left blank.

jpd_percent = cwmat.freq;
jpd_percent.values = 100 * cwmat.freq.values;
jpd_percent.values(jpd_percent.values == 0) = NaN;
figure('Position', [100 100 1000 800]);
ax2 = mhkit_plot_matrix(jpd_percent, 'xlabel', 'Te [s]', 'ylabel', 'Hm0 [m]', ...
    'zlabel', 'Occurrence [%]', 'trim_to_data', true, 'value_format', '%.2f %%', 'font_size', 8);
title(ax2, sprintf('Joint Probability Distribution | %s', title_suffix));
%%
% The underlying |mhkit_plot_matrix| function accepts any matrix struct and
% custom axis labels. Here the mean power matrix is plotted in kW.

avg_power_mat_kW = avg_power_mat;
avg_power_mat_kW.values = avg_power_mat.values / 1000;
figure('Position', [100 100 1000 800]);
ax3 = mhkit_plot_matrix(avg_power_mat_kW, 'xlabel', 'Te [s]', 'ylabel', 'Hm0 [m]', ...
    'zlabel', 'Mean Power [kW]', 'trim_to_data', true, 'value_format', '%.2f kW', 'font_size', 8);
title(ax3, sprintf('Mean Power Matrix | %s', title_suffix));
