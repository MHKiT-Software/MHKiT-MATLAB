%% MHKiT WPTO Hindcast Example
% This example demonstrates how to load and analyze data from the WPTO hindcast
% dataset hosted on AWS.
%
% <https://registry.opendata.aws/wpto-pds-us-wave/ WPTO Hindcast Dataset>

%% Dataset Description
% The dataset consists of two main components:
%
% *High-spatial-resolution dataset ('3-hour'):*
%
% * Covers U.S. Exclusive Economic Zone (EEZ) along West coast and Hawaii
% * Uses unstructured grid with ~200m resolution in shallow water
% * 3-hour time step resolution
% * Spans 32 years (01/01/1979 - 12/31/2010)
%
% *Virtual buoy dataset ('1-hour'):*
%
% * Available at specific locations within the spatial domain
% * Spans same 32-year period
% * 1-hour time resolution

%% Available Variables
% *3-hour Dataset Variables:*
%
% Variables are indexed by latitude, longitude, and time:
%
% * |significant_wave_height|: significant wave height, $H_{m_0}$ [m]
% * |energy_period|: energy period, $T_e$ [s]
% * |peak_period|: peak period, $T_p$ [s]
% * |mean_absolute_period|: mean absolute period, $T_m$ [s]
% * |mean_zero-crossing_period|: mean zero-crossing period, $T_z$ [s]
% * |omni-directional_wave_power|: omnidirectional wave power, $J$ [W/m]
% * |spectral_width|: spectral width, $\epsilon_0$ [-]
% * |directionality_coefficient|: directionality coefficient, $d$ [-]
% * |maximum_energy_direction|: direction of maximum wave energy, $\theta_{J_{max}}$ [deg]
% * |mean_wave_direction|: mean wave direction, $\theta_m$ [deg]
% * |water_depth|: water depth, $h$ [m]
%
% *1-hour Dataset Variables:*
%
% Includes all variables from 3-hour dataset plus:
%
% * |directional_wave_spectrum|: directional wave spectrum, $S(f, \theta)$
% * |frequency_bin_edges|: frequency bin edges, $f$ [Hz]

%% Data Access Configuration
%
% To access the WPTO hindcast data:
%
% # Obtain API key from <https://developer.nlr.gov/signup/>

%% Example 1: Request Single Location Data
%
% Location: <https://www.energy.gov/eere/water/pacwave-offshore-wave-energy-test-site PacWave South>
%
% Request 3-hour significant wave height, $H_{m_0}$ [m], data for 1995 at PacWave South
% This example demonstrates basic data retrieval for a single parameter at one location

% Set parameters for data request
data_type = '3-hour';
year = 1995;
lat_lon = [44.624076, -124.280097]; % PacWave South
parameter = ["significant_wave_height"];  % Requested parameter
api_key = '3K3JQbjZmWctY0xmIfSYvYgtIcM3CN0cb1Y2w9bf';  % Demo API key (rate-limited)

% Request data for the specified location and time period
wave_data = request_wpto(data_type, parameter, lat_lon, year, api_key);

%%
% Plot the significant wave height, $H_{m_0}$ [m], time series

figure('Position', [100, 100, 1000, 400]);
plot(wave_data.time, wave_data.significant_wave_height, 'LineWidth', 1.5);
title(['Significant Wave Height, H_{m_0}, at ' num2str(wave_data.metadata.latitude) '°N, ' ...
       num2str(-1 * wave_data.metadata.longitude) '°W']);
xlabel('Time [UTC]');
ylabel('Significant Wave Height, H_{m_0} [m]');
grid on;
xtickformat('MMM-yy');

%% Example 2: Request Multiple Locations and Parameters
%
% Locations:
%
% * <https://www.energy.gov/eere/water/pacwave-offshore-wave-energy-test-site PacWave South>
% * <https://www.energy.gov/eere/water/pacwave-offshore-wave-energy-test-site PacWave North>
%
% Request 3-hour energy period, $T_e$ [s], and significant wave height, $H_{m_0}$ [m], at two locations
% This example shows how to handle multiple parameters and locations simultaneously

% Define multiple parameters and locations
parameter = ["energy_period", "significant_wave_height"];
lat_lon = [44.624076, -124.280097; % PacWave South
           43.489171, -125.152137]; % PacWave North

% Request data for both locations
wave_measurements = request_wpto(data_type, parameter, lat_lon, year, api_key);

%%
% Create subplots for energy period, $T_e$, and significant wave height, $H_{m_0}$, at both locations

figure('Position', [100, 100, 1200, 800]);

% Plot Energy Period
subplot(2,1,1);
plot(wave_measurements(1).time, wave_measurements(1).energy_period, '-', 'LineWidth', 1.5);
hold on;
plot(wave_measurements(2).time, wave_measurements(2).energy_period, '--', 'LineWidth', 1.5);
hold off;
title('Energy Period, T_e, at Two Locations');
xlabel('Time [UTC]');
ylabel('Energy Period, T_e [s]');
grid on;
legend(['Location 1 (' num2str(wave_measurements(1).metadata.latitude) '°N, ' ...
        num2str(-1 * wave_measurements(1).metadata.longitude) '°W)'], ...
       ['Location 2 (' num2str(wave_measurements(2).metadata.latitude) '°N, ' ...
        num2str(-1 * wave_measurements(2).metadata.longitude) '°W)']);
xtickformat('MMM-yy');

% Plot Significant Wave Height
subplot(2,1,2);
plot(wave_measurements(1).time, wave_measurements(1).significant_wave_height, '-', 'LineWidth', 1.5);
hold on;
plot(wave_measurements(2).time, wave_measurements(2).significant_wave_height, '--', 'LineWidth', 1.5);
hold off;
title('Significant Wave Height, H_{m_0}, at Two Locations');
xlabel('Time [UTC]');
ylabel('Significant Wave Height, H_{m_0} [m]');
grid on;
legend(['Location 1 (' num2str(wave_measurements(1).metadata.latitude) '°N, ' ...
        num2str(-1 * wave_measurements(1).metadata.longitude) '°W)'], ...
       ['Location 2 (' num2str(wave_measurements(2).metadata.latitude) '°N, ' ...
        num2str(-1 * wave_measurements(2).metadata.longitude) '°W)']);
xtickformat('MMM-yy');

%% Example 3: Request Peak Period and Wave Power
%
% Location: <https://www.energy.gov/eere/water/pacwave-offshore-wave-energy-test-site PacWave South>
%
% Request 3-hour peak period, $T_p$ [s], and omnidirectional wave power, $J$ [W/m], data
% This example shows the wave resource quantities used to characterize a site

% Set parameters for the data request
data_type = '3-hour';
parameter = ["peak_period", "omni-directional_wave_power"];
lat_lon = [44.624076, -124.280097]; % PacWave South

% Request data
resource_data = request_wpto(data_type, parameter, lat_lon, year, api_key);

%%
% Plot the peak period, $T_p$ [s], and omnidirectional wave power, $J$. The hindcast
% stores $J$ in [W/m], it is plotted in [kW/m]

% Peak Period Plot
figure('Position', [100, 100, 1000, 400]);
plot(resource_data.time, resource_data.peak_period, 'LineWidth', 1.5);
title(['Peak Period, T_p, at ' num2str(resource_data.metadata.latitude) '°N, ' ...
       num2str(-1 * resource_data.metadata.longitude) '°W']);
xlabel('Time [UTC]');
ylabel('Peak Period, T_p [s]');
grid on;
xtickformat('MMM-yy');

% Omnidirectional Wave Power Plot
figure('Position', [100, 100, 1000, 400]);
plot(resource_data.time, resource_data.omni_directional_wave_power / 1000, 'LineWidth', 1.5);
title(['Omnidirectional Wave Power, J, at ' num2str(resource_data.metadata.latitude) '°N, ' ...
       num2str(-1 * resource_data.metadata.longitude) '°W']);
xlabel('Time [UTC]');
ylabel('Omnidirectional Wave Power, J [kW/m]');
grid on;
xtickformat('MMM-yy');
