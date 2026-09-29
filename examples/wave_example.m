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
%% Joint Probability Distribution
% The joint probability distribution (JPD) is the fraction of records that
% fall in each Hm0 and Te bin, and it depends only on the wave resource.
% |plot_wave_joint_probability_distribution| bins the sea states on the
% IEC/TS 62600-100 grid (0.5 m Hm0 and 1 s Te bins with edges on 0, 0.5, 1,
% ... m and 0, 1, 2, ... s), plots the JPD as a percentage with the empty bins
% left blank, and returns the JPD matrix struct.

figure('Position', [100 100 1000 800]);
[~, jpd] = plot_wave_joint_probability_distribution(Hm0, Te, ...
    'title', sprintf('Joint Probability Distribution | %s', title_suffix));

% The same bins are used for the capture width and power matrices below
Hm0_bins = jpd.y_bins;
Te_bins = jpd.x_bins;
%% Model Device Power Using DOE Reference Models
% In a real application, the user would provide power values measured from a
% WEC. This example instead looks up each sea state in the power matrix of a
% device from the <https://openei.org/wiki/PRIMRE/Signature_Projects/Reference_Model
% Reference Model Project>, a U.S. Department of Energy effort that developed
% open-source marine energy point designs to benchmark technology performance
% and costs. MHKiT ships the matrices for Reference Model 3 (wave point
% absorber), 5 (oscillating surge flap), and 6 (oscillating water column) from
% the <https://github.com/NatLabRockies/SAM National Laboratory of the Rockies
% (NLR) System Advisor Model (SAM)>; change |device| to switch. The devices are
% described in Neary et al. (2014), <https://doi.org/10.2172/1159756 Methodology
% for Design and Economic Analysis of Marine Energy Conversion (MEC) Technologies>,
% SAND2014-9040, Sandia National Laboratories.

device = "RM3";  % "RM3", "RM5", or "RM6"
reference_model = reference_model_power_matrix(device);
disp(reference_model.description)
disp(reference_model.data_source)

reference_model_kW = reference_model;
reference_model_kW.values = reference_model.values / 1000;
figure('Position', [100 100 1000 800]);
ax_device = mhkit_plot_matrix(reference_model_kW, ...
    'xlabel', 'Energy Period, T_e [sec]', 'ylabel', 'Significant Wave Height, H_{m0} [m]', ...
    'zlabel', 'Power [kW]', 'trim_to_data', true, 'value_format', '%.0f', 'font_size', 8);
title(ax_device, sprintf('%s Power Matrix | Neary et al. (2014)', device));

% Interpolate between bin centers; sea states outside the matrix produce no power
Power = interp2(reference_model.x_bins, reference_model.y_bins, reference_model.values, ...
    Te, Hm0, 'linear', 0);  % [W]

figure('Position', [100 100 900 350]);
plot(time, Power / 1000); ylabel('Modeled Power [kW]'); xlabel('Time [UTC]'); grid on
title(sprintf('Modeled %s Power | %s', device, title_suffix))
%% Capture Width
% The following operations create capture width matrices, as specified by the
% IEC/TS 62600-100, on the same Hm0 and Te bins as the JPD. But first, we need
% to calculate capture width. Keep in mind that the power comes from a modeled
% power matrix, not measurements, so the capture width varies smoothly with
% the sea state.

% calculating capture width with power and wave flux in vectors
CW = capture_width(Power, J);

% Calculate the necessary capture width matrices for each statistic based
% on IEC/TS 62600-100
cwmat.mean = capture_width_matrix(Hm0, Te, CW, "mean", Hm0_bins, Te_bins);
cwmat.std = capture_width_matrix(Hm0, Te, CW, "std", Hm0_bins, Te_bins);
cwmat.count = capture_width_matrix(Hm0, Te, CW, "count", Hm0_bins, Te_bins);
cwmat.min = capture_width_matrix(Hm0, Te, CW, "min", Hm0_bins, Te_bins);
cwmat.max = capture_width_matrix(Hm0, Te, CW, "max", Hm0_bins, Te_bins);
%%
% The graphics function |plot_matrix| visualizes a matrix struct. Each bin is
% drawn between its edges with the tick marks at the bin edges, and
% |trim_to_data| limits the axes to the bins that hold data plus one empty
% bin around them.

figure('Position', [100 100 1000 800]);
ax_cw = plot_matrix(cwmat.mean, "Capture Width", "zlabel", "Capture Width [m]", ...
    "trim_to_data", true, "font_size", 8);
title(ax_cw, sprintf('Modeled %s Mean Capture Width Matrix | %s', device, title_suffix));
%% Power Matrix
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
%%
% The underlying |mhkit_plot_matrix| function accepts any matrix struct and
% custom axis labels. Here the mean power matrix is plotted in kW.

avg_power_mat_kW = avg_power_mat;
avg_power_mat_kW.values = avg_power_mat.values / 1000;
figure('Position', [100 100 1000 800]);
ax_power = mhkit_plot_matrix(avg_power_mat_kW, 'xlabel', 'Energy Period, T_e [sec]', 'ylabel', 'Significant Wave Height, H_{m0} [m]', ...
    'zlabel', 'Mean Power [kW]', 'trim_to_data', true, 'value_format', '%.1f kW', 'font_size', 8);
title(ax_power, sprintf('Modeled %s Mean Power Matrix | %s', device, title_suffix));
%% Energy Production
% The modeled power series gives the energy production of the device at this
% site, first month by month and then as the mean annual energy production
% (MAEP) defined by IEC/TS 62600-100.
%% Monthly Energy Production
% Each record is a 30 minute average, so the energy per record is the power
% times that duration. Summing the records by month shows how the energy
% production is distributed through the year.

record_hours = hours(median(diff(time)));
energy_MWh = Power * record_hours / 1e6;  % [MWh] per record
monthly_MWh = accumarray(month(time), energy_MWh, [12 1]);

figure('Position', [100 100 900 400]);
bar(1:12, monthly_MWh);
xticks(1:12); xticklabels(month(datetime(2025, 1:12, 1), 'shortname'));
ylabel('Energy [MWh]'); grid on
yline(mean(monthly_MWh), '--', 'Monthly mean', 'LabelHorizontalAlignment', 'left');
title(sprintf('Modeled %s Monthly Energy Production | %s', device, title_suffix))
%% Mean Annual Energy Production
% There are two ways to calculate mean annual energy production (MAEP). One
% is from capture width and wave energy flux matrices, the other is from time
% series data, as shown below.

% Calculate maep from timeseries
maep_timeseries = mean_annual_energy_production_timeseries(CW, J);  % [W*h]
% Calculate maep from matrix
maep_matrix = mean_annual_energy_production_matrix(cwmat.mean, jmat, jpd);  % [W*h]

% MHKiT returns energy in W*h; report in MWh
fprintf('Modeled %s MAEP from time series: %.0f MWh\n', device, maep_timeseries / 1e6);
fprintf('Modeled %s MAEP from matrices:    %.0f MWh\n', device, maep_matrix / 1e6);
