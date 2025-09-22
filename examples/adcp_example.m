%% Analyzing ADCP Data with DOLfYN
%
% The following example walks through a straightforward workflow for analyzing
% Acoustic Doppler Current Profiler (ADCP) data utilizing MHKiT-MATLAB.
%
% MHKiT-MATLAB has developed native MATLAB code based on the Doppler
% Oceanographic Library for pYthoN (DOLfYN) in MHKiT-Python as a module to
% facilitate ADCP and Acoustic Doppler Velocimetry (ADV) data processing.
%
% Here is a standard workflow for ADCP data analysis:
%
%
% *1. Import Data*
%
% *2. Review, QC, and Prepare the Raw Data:*
%
% # Calculate or verify the correctness of depth bin locations
% # Discard data recorded above the water surface or below the seafloor
% # Assess the quality of velocity, beam amplitude, and/or beam correlation data
% # Rotate Data Coordinate System
%
% *3. Data Averaging:*
%
% * If not already executed within the instrument, average the data into time bins of a predetermined duration, typically between 5 and 10 minutes
%
% *4. Calculation of Speed and Direction*
%
% *5. Visualizations*
%
% *6. Saving and Loading DOLfYN datasets*
%
% *7. Turbulence Statistics:*
%
% * Turbulence Intensity (TI)
% * Power Spectral Densities
% * Instrument Noise
% * Turbulent Kinetic Energy (TKE) Dissipation Rate
% * Noise-corrected TI
% * TKE Components
% * TKE Production
% * TKE Balance

%% 1. Importing Raw Instrument Data
%
% One of DOLfYN's key features is its ability to directly import raw data from
% an Acoustic Doppler Current Profiler (ADCP) right after it has been
% transferred. This example uses a Nortek Signature1000 ADCP, with
% the data stored in files with an '.ad2cp' extension. This specific dataset
% represents several hours of velocity data, captured at 1 Hz by an ADCP
% mounted on a bottom lander within a tidal inlet. The list of instruments
% compatible with DOLfYN can be found in the
% <https://mhkit-software.github.io/MHKiT/mhkit-python/api.dolfyn.html MHKiT DOLfYN documentation>.
%
% To start, this section uses the |dolfyn_read| function to import the raw ADCP data into a MATLAB struct.
% The |dolfyn_read| function processes the input raw file and converts the data
% into a MATLAB struct, typically called a dataset, or |ds| for short. This
% dataset includes several groups of variables:
%
% # *Velocity*: Recorded in the coordinate system saved by the instrument (beam, XYZ, ENU)
% # *Beam Data*: Includes amplitude and correlation data
% # *Instrumental & Environmental Measurements*: Captures the instrument's bearing and environmental conditions
% # *Orientation Matrices*: Used by DOLfYN for rotating through different coordinate frames.
%
% Note: Before reading ADCP data with |dolfyn_read|, ensure your files are split
% into manageable chunks. DOLfYN performs best with files under 1GB. The file
% size depends primarily on your sampling frequency and deployment duration.
% Most ADCP manufacturers provide software tools or deployment settings to
% control output file size. Consult your instrument's documentation for
% specific guidance on file size configuration.
%
% The |dolfyn_read| function parses the raw ADCP file into MATLAB.

ds = dolfyn_read('./examples/data/dolfyn/Sig1000_tidal.ad2cp');

%%
%
% The |dolfyn_read| function automatically detects the instrument type and reads
% the specified ensemble data (by default all data) into a MATLAB struct.
% ADCP datasets contain multiple measurement types, coordinate systems, and data dimensions
% organized in a hierarchical structure. When working with raw ADCP data, it is recommended to
% examine the dataset organization before proceeding with analysis. Understanding how the ADCP data is stored
% will add value to both the analysis and interpretation of the data.
%
% The |ds| struct contains many fields. The following key fields define important aspects of the dataset:

%%
% |attrs| contains global attributes that provide metadata about the dataset,
% including instrument specifications, deployment configuration, and processing history.

ds.attrs

%%
% |coords| contains the coordinate arrays that define the dimensional structure of the dataset. These coordinates specify how multidimensional data arrays are indexed and organized.

ds.coords

%%
% |vel| contains velocity measurements organized by dimensions and coordinates,
% with the actual measurement data stored in the |data| field.

ds.vel

%%
% |ds.vel.dims| shows the dimension labels of the velocity data that correspond
% to the actual shape of the underlying data array.

ds.vel.dims
size(ds.vel.data)

%%
% The output shows that the velocity data contains 55,000 time samples,
% 28 range bins, and 4 directional components (typically the 4 beam velocities).
% All MHKiT |dolfyn_| functions are designed to work with these dimensions and
% perform necessary data reshaping internally.

%% Using |dolfyn_plot| to visualize ADCP Data:
%
% The |dolfyn_plot| function is designed to visualize multidimensional ADCP data directly from DOLfYN datasets.
% Pass in the dataset, the dimension name and index.
% |dolfyn_plot| uses data found within the dataset attributes to determine
% the plot type and labels.
%
% The following visualization parameters define standard figure dimensions
% for consistent plot scaling. These viz_width and viz_height variables
% are passed to |dolfyn_plot| to maintain uniform visualization sizing.

viz_width = 1800;
viz_height = 300;

%%
% |dolfyn_plot| is a versatile visualization function for DOLfYN datasets. This example uses it
% to create histograms of velocity measurements across all depth bins. The function
% accepts the following key arguments:
%
% * ds: the dataset structure
% * 'vel': specifies the variable in |ds| to plot
% * 'dim', specifies the dimension to plot, 1:3: plots all three velocity components (typically east, north, up)
% * 'width'/'height': control figure dimensions
% * 'kind', 'hist': specifies histogram plots instead of the default time series
%
% The resulting histograms help identify potential outliers and assess the overall
% distribution of velocity measurements.

dolfyn_plot(ds, 'vel', 'dim', 1:3, 'width', viz_width, 'height', viz_height, 'kind', 'hist');

%%
% Plot the data with a logarithmic y scale to identify outliers, particularly in the extremes of the dataset.

dolfyn_plot(ds, 'vel', 'dim', 1:3, 'width', viz_width, 'height', viz_height, 'kind', 'hist', 'scale', 'logy');

%%
% For this dataset |./examples/data/dolfyn/Signature1000_tidal.ad2cp| the
% histogram plots show that the velocity measurements are symmetrically
% distributed and centered around 0 m/s, with no significant spikes or
% discontinuities at the extremes. This indicates reasonable data quality
% and the analysis can proceed without removing outliers. Unusual peaks at extreme
% values or highly asymmetric distributions would require further investigation.
%
% Note: While visual inspection provides an initial assessment, a complete
% quality control process would typically include additional statistical
% tests and instrument-specific checks. This example proceeds
% with the original data.

%%
% To ensure consistent scaling across velocity plots, this code calculates and defines the colorbar limits
% The velocity colorbar limits are determined using percentile-based thresholds
% to reduce the impact of outliers while maintaining symmetry around zero.
% First, the 1st and 99th percentiles of the velocity data are calculated,
% capturing the central 98% of the data distribution. Then, the larger absolute
% value between these percentiles is used to set symmetric positive and negative
% limits. The final values are rounded to whole numbers for a cleaner colorbar
% scale. This approach ensures the colorbar:
%
% # Excludes extreme outliers that might compress the color scale
% # Maintains symmetry around zero for proper visualization of bidirectional flow
% # Preserves the true scale of the dominant velocity magnitudes
% # Uses whole numbers for easier readability

floor_percentile = 1;
ceiling_percentile = 99;

% Calculate 1st and 99th percentiles of the flattened velocity array
pct_min = prctile(ds.vel.data(:), floor_percentile);
pct_max = prctile(ds.vel.data(:), ceiling_percentile);

% Determine the larger magnitude between percentile bounds
max_abs = max(abs([pct_min, pct_max]));

% Set symmetric colorbar limits using the larger magnitude, rounded to whole numbers
vel_cbar_min = -ceil(max_abs);
vel_cbar_max = ceil(max_abs);

%%
% The following creates a visualization of the ADCP velocity data.
% dolfyn_plot will generate three panels showing:
%
% # East velocity component (dim=1): positive values indicate eastward flow
% # North velocity component (dim=2): positive values indicate northward flow
% # Up velocity component (dim=3): positive values indicate upward flow
%
% Each panel shows:
%
% * X-axis: Time
% * Y-axis: Distance from ADCP (range)
% * Colors: Velocity magnitude (in m/s)
% * Blues/Negatives: flow in negative direction
% * Reds/Positives: flow in positive direction
%
% The previously calculated |vel_cbar_min| and |vel_cbar_max| ensure
% color scales are consistent and centered around zero across all plots.

dolfyn_plot(ds, 'vel', 'dim', 1:3, ...
    'width', viz_width, 'height', viz_height, ...
    'cbar_min', vel_cbar_min, 'cbar_max', vel_cbar_max);

%%
% Velocity dimensions are located in |ds.vel.dims|

ds.vel.dims

%%
% The specifications for each coordinate are located in |ds.coords.(dim)|

ds.coords.dir

%% 2. Initial Steps for Data Quality Control (QC)
%
%
% *2.1: Set the Deployment Height*
%
% When using Nortek instruments, the deployment software does not factor in the
% deployment height. The deployment height represents the position of the
% Acoustic Doppler Current Profiler (ADCP) within the water column.
%
% In this context, the center of the first depth bin is situated at a distance
% that is the sum of three elements:
%
% # Deployment height (the ADCP's position in the water column)
% # Blanking distance (the minimum distance from the ADCP to the first measurement point)
% # Cell size (the vertical distance of each measurement bin in the water column)
%
% To ensure accurate readings, it is critical to calibrate the 'range'
% coordinate to make '0' correspond to the seafloor. This calibration can be
% achieved using the |set_range_offset| function. This function is also useful
% when working with a down-facing instrument as it helps account for
% the depth below the water surface.
%
% For Teledyne RDI ADCPs, the TRDI deployment software will prompt
% the user to specify the deployment height/depth during setup. If there's a need
% for calibration post-deployment, the |set_range_offset| function can be
% utilized in the same way as described above.

%%
% The ADCP transducers were measured to be 0.6 meters from the feet of the lander

ds = set_range_offset(ds, 0.6);

%%
% Range values can be checked by calculating min, max, and spacing values.

% Create a table summarizing the range characteristics
range_stats = table({"Maximum"; "Minimum"; "Mean Spacing"; "Max Spacing"; "Min Spacing"}, ...
    round([max(ds.coords.range); ...
          min(ds.coords.range); ...
          mean(diff(ds.coords.range)); ...
          max(diff(ds.coords.range)); ...
          min(diff(ds.coords.range))] * 100) / 100, ...
    'VariableNames', {'Range Metric', 'Value'});

format short;
range_stats

%%
% *2.2. Discard Data Above Surface Level*
%
% To reduce computational load, this section excludes all data at or above the water
% surface level. Since the instrument was oriented upwards, this approach utilizes the
% pressure sensor data along with the function |water_depth_from_pressure|. However,
% this approach necessitates that the pressure sensor was calibrated or
% 'zeroed' prior to deployment. If the instrument is facing downwards or
% doesn't include pressure data, the function |water_depth_from_pressure| can be used to
% detect the seabed or water surface.
%
% It's important to note that Acoustic Doppler Current Profilers (ADCPs) do not
% measure water salinity, so the user will need to supply this information to
% the function. The dataset returned by this function includes an additional
% variable, "depth". If |water_depth_from_pressure| is invoked after
% |set_range_offset|, "depth" represents the distance from the water
% surface to the seafloor. Otherwise, it indicates the distance to the
% ADCP pressure sensor. The |water_depth_from_pressure| function allows
% the user to specify a salinity value in Practical Salinity
% Units (PSU).

water_salinity_psu = 31;
ds = water_depth_from_pressure(ds, 'salinity', water_salinity_psu);

%%
% Visualizing Depth

dolfyn_plot(ds, 'depth', 'width', 800, 'height', 400);
title('Distance from ADCP Pressure Sensor (Including Offset) to Water Surface [m]');

%%
% After determining the depth/altitude, the |remove_surface_interference| function can be used to
% discard data in depth bins near or above the actual water surface. This function
% calculates a range limit based on the beam angle and cell size, where surface
% interference is expected to occur at distances greater than |range * cos(beam angle)
% - cell size|. The beam angle accounts for the acoustic spread of the ADCP signal
% (typically 20-25 degrees). To exclude the found surface interference data DOLfYN sets data above this threshold to NaN.

ds = remove_surface_interference(ds);

%%
% Visualizing removal of surface interference

dolfyn_plot(ds, 'vel', 'dim', 1:3, ...
    'width', viz_width, 'height', viz_height, ...
    'cbar_min', vel_cbar_min, 'cbar_max', vel_cbar_max);

%%
% *Correlation*
%
% Reviewing data from the other beams is recommended. A significant
% portion of this data is of high quality. To avoid discarding valuable data
% with lower correlations, which could be due to natural variations, this section uses
% the |correlation_filter| function to assign a value of NaN (not a number)
% to velocity values corresponding to correlations below 50%.
%
% However, it's important to note that the correlation threshold is dependent
% on the specifics of the deployment environment and the instrument used. It's
% not unusual to set a threshold as low as 30%, or even to forgo the use of
% this function entirely.

ds = correlation_filter(ds, 'thresh', 50);

%%
% Visualizing Correlation Data

dolfyn_plot(ds, 'corr', 'dim', 1:4, 'width', viz_width, 'height', viz_height * 2);

%%
% *2.4 Rotate Data Coordinate System*
%
% After cleaning the data, the next step is to rotate the velocity data into accurate East, North, Up (ENU) coordinates.
%
% ADCPs utilize an internal compass or magnetometer to determine magnetic ENU
% directions. You can use the set_declination function to adjust the velocity
% data according to the magnetic declination specific to your geographical
% coordinates. This declination can be looked up online for specific
% coordinates.
%
% Instruments save vector data in the coordinate system defined in the
% deployment configuration file. To make this data meaningful, it must be
% transformed through various coordinate systems
% ("beam"<->"inst"<->"earth"<->"principal"). This transformation is
% accomplished using the |rotate2| function. If the "earth" (ENU) coordinate
% system is specified, DOLfYN will automatically rotate the dataset through
% the required coordinate systems to reach the "earth" coordinates.
%
% In this case, since the ADCP data is already in the "earth" coordinate
% system, the |rotate2| function will return the input dataset without
% modifications. The |set_declination| function will work no matter the
% coordinate system.

ds = set_declination(ds, 15.8);
ds = rotate2(ds, 'earth');

%%
% For analysis that requires changing the frame of reference, DOLfYN provides
% the |rotate2| function to transform data between coordinate systems
% (beam, inst, earth, principal). To rotate into the principal frame
% (streamwise, cross-stream, vertical), first calculate the depth-averaged
% principal flow heading and add it to the dataset attributes.

ds.attrs.principal_heading = calc_principal_heading(squeeze(ds.vel.data), true);
disp(ds.attrs.principal_heading);
ds_streamwise = rotate2(ds, 'principal');

%%
% Visualize Streamwise Velocity

%%
% First verify the transformation by visually inspecting the histogram and checking for outliers

dolfyn_plot(ds_streamwise, 'vel', 'dim', 1:3, 'width', viz_width, 'height', viz_height, 'kind', 'hist');

%%
% Next plot streamwise and cross-stream velocity

dolfyn_plot(ds_streamwise, 'vel', 'dim', 1:2, ...
    'width', viz_width, 'height', viz_height, ...
    'cbar_min', vel_cbar_min, 'cbar_max', vel_cbar_max);

%% 3. Average the Data
%
% As this deployment was configured in "burst mode", a standard step in the
% analysis process is to average the velocity data into time bins.
%
% However, if the instrument was set up in an "averaging mode" (where a
% specific profile and/or average interval was set, for instance, averaging
% 5 minutes of data every 30 minutes), this step would have been
% performed within the ADCP during deployment and can thus be skipped.
%
% The function |average_by_dimension| is used to average the data by dimension.
% To average the data into time bins (also known as ensembles), the sampling rate and averaging period must be defined.
% DOLfYN stores the sampling frequency [Hz] in |.attrs.fs|. The averaging period should be specified in seconds.
% These values specify the number of samples to average.

sampling_frequency_hz = ds.attrs.fs;
averaging_period_seconds = 300; % 5 minutes

averaging_samples = int32(round(averaging_period_seconds * sampling_frequency_hz));

%%
% Once the averaging samples count has been calculated, the
% dataset can be averaged by time by passing the original dataset, the number of samples to
% average, and the dimension ('time').
%
% Important note about sample size:
% When the total number of samples is not evenly divisible by your averaging period,
% DOLfYN will discard the remaining samples at the end of the dataset. For example,
% with:
%
% * Total samples: 55000
% * Averaging period: 300 samples
% * 55000 ÷ 300 = 183.33 (not a whole number)
% * Last 100 samples will be discarded
%
% To avoid losing data, you can:
%
% # Combine multiple datasets together
% # Adjust your averaging period to evenly divide into your total samples
% # Trim your dataset to a length that's divisible by your averaging period
% # Accept the loss of the partial bin at the end if it's not critical
%
% In this example, losing 100 samples (100/55000 = 0.18% of data) is acceptable
% for this analysis.

ds_avg = average_by_dimension(ds, averaging_samples, 'time');

%%
% Plotting the Summary Data

dolfyn_plot(ds_avg, 'vel', 'dim', 1:3, 'width', viz_width, 'height', viz_height, 'kind', 'hist');

dolfyn_plot(ds_avg, 'vel', 'dim', 1:3, ...
    'width', viz_width, 'height', viz_height, ...
    'cbar_min', vel_cbar_min, 'cbar_max', vel_cbar_max);

%% 4. Calculate Current Speed and Direction
%
% The |calculate_horizontal_speed_and_direction| function processes the east and
% north velocity components to compute two value added fields: current speed (U_mag)
% and direction (U_dir). The speed is calculated as the magnitude of the horizontal
% velocity components and stored in ds_avg.U_mag, using the same units as the
% input velocities (typically m/s).
%
% Direction (U_dir) handling depends on the coordinate system (ds.coord_sys):
%
% * For "earth" coordinates: Direction is measured in degrees clockwise from true North [0-360°]
% * For other coordinates (e.g., "principal", "inst", "beam"): Direction is measured in degrees clockwise from East/streamwise component [0-360°]
%
% For example in "earth" coordinates:
%
% * U_dir = 90° means the current is flowing towards the east
% * U_dir = 180° means the current is flowing towards the south
% * U_dir = 270° means the current is flowing towards the west
%
% See 'help calculate_horizontal_speed_and_direction' for complete details on
% coordinate systems and direction calculations.

ds_avg = calculate_horizontal_speed_and_direction(ds_avg);

ds_avg.U_mag
ds_avg.U_dir

%% 5. Visualize Current Speed and Direction
%
% The following creates two visualizations to illustrate the flow patterns:
%
% *1. Current Speed Plot:*
%
% * Shows the magnitude of water movement at each depth and time
% * Uses a blue colormap where darker colors indicate faster flow
% * Includes a blue line showing the water surface elevation
% * X-axis: Time of day
% * Y-axis: Distance from seafloor (altitude)
%
% *2. Current Direction Plot:*
%
% * Shows where the water is flowing TO at each depth and time
% * Uses a circular colormap to represent the full 360° of possible directions
% * Includes the same water surface elevation line
% * Same axes as the speed plot for direct comparison
%
% Both plots use the time-averaged data (ds_avg), where each cell represents a 5-minute average.

% Create color map for both visualizations
blues = cmocean('ice-', 256);
surface_elevation_blue = blues(128,:);

%%
% Current Speed Visualization

fig = figure;
viz_width = 800;
fig.Position = [100 100 viz_width viz_height];

ax = axes('Position', [0.14 0.14 0.8 0.74]);

time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
range = ds_avg.coords.range;
speed = squeeze(ds_avg.U_mag.data);
speed = speed';
depth = ds_avg.depth.data;
y_min = 0;
% Calculate y_max so all plots have the same y-axis limits
y_max = int32(max(ds.depth.data) + 1);

[T, R] = meshgrid(time, range);
pcolor(T, R, speed);
shading flat;

colormap(blues);

% Ensure that the box and ticks are on top of the visualization
set(gca, 'Layer', 'top', 'Box', 'on')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
set(gca, 'TickDir', 'in')
set(gca, 'LineWidth', 1)

hold on;

% Add a line denoting water surface depth using the darkest blue from the colormap
surface_elevation = plot(time, depth, 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'southeast');

xlabel('Time');
ylabel('Altitude [m]');
ylim([y_min y_max]);

ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = sprintf('%s [%s]', ds_avg.U_mag.long_name, ds_avg.U_mag.units);
title('Water Speed [m/s]');

hold off;

%%
% Current "To" Direction (Degrees CW from True North) Visualization

fig = figure;
fig.Position = [100 100 viz_width viz_height];

ax = axes('Position', [0.14 0.14 0.8 0.74]);

direction = squeeze(ds_avg.U_dir.data);
direction = direction';

pcolor(T, R, direction);
shading flat;
colormap(cmocean('phase', 256));

set(gca, 'Layer', 'top', 'Box', 'on')
set(gca, 'XGrid', 'off', 'YGrid', 'off')
set(gca, 'TickDir', 'in')
set(gca, 'LineWidth', 1)

hold on;

% Plot water surface depth
surface_elevation = plot(time, depth, 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'southeast');

xlabel('Time');
ylabel('Altitude [m]');
ylim([y_min y_max]);

ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = "Horizontal Vel Dir [deg CW from true N]";

title('Horizontal Velocity "To" Direction [deg]');

hold off;

%% 6. Save Data
%
% Datasets can be saved using the dolfyn |write_netcdf| function and be read
% using the dolfyn |read_netcdf| function.

% Uncomment these lines to save and load to your current working directory
% write_netcdf(ds, 'your_data.nc');
% ds_saved = read_netcdf('your_data.nc');

%% 7. Turbulence Statistics
%
% The next section of this example will run through the turbulence analysis 
% of the data presented here. There was no intention of measuring turbulence 
% in the deployment that collected this data, so results depicted here are 
% not the highest quality. The quality of turbulence measurements from an 
% ADCP depend heavily on the quality of the deployment setup and data 
% collection, particularly instrument frequency, sampling frequency and 
% depth bin size.
%
% Read more on proper ADCP setup for turbulence measurements in: Thomson, Jim,
% et al. "Measurements of turbulence at two tidal energy sites in Puget
% Sound, WA." IEEE Journal of Oceanic Engineering 37.3 (2012): 363-374.
% <https://doi.org/10.1109/JOE.2012.2191656>
%
% Most functions related to turbulence statistics in MHKiT-MATLAB dolfyn module
% have the papers they originate from referenced in their docstrings.

%% 7.1 Turbulence Intensity
% For most users, turbulence intensity (TI), the ratio of the ensemble
% standard deviation to ensemble flow speed given as a percent, is the primary
% metric needed. In MHKiT, this can be simply calculated from ensemble-averaged 
% data, but be aware that this will be a conservative estimate. Another 
% function, |calculate_turbulence_intensity|, is capable of subtracting 
% instrument noise from this parameter and is discussed below. The 
% noise-subtracted TI is more accurate and typically 1-2% lower than the 
% non-noise-subtracted estimation.

% Basic turbulence intensity calculation
ds_avg = calculate_turbulence_intensity(ds_avg, 'field_name', 'TI');

% Create visualization
fig = figure;
fig.Position = [100 100 800 400];

ti_data = ds_avg.TI.data';  % Transpose for proper orientation
time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
range = ds_avg.coords.range;
depth = ds_avg.depth.data;

[T, R] = meshgrid(time, range);
pcolor(T, R, ti_data);
shading flat;
colormap(cmocean('amp', 256));  % Red colormap similar to Python

hold on;

% Plot water surface depth
surface_elevation = plot(time, depth, 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'southeast');

xlabel('Time');
ylabel('range');
ylim([0, max(depth) + 1]);

ax = gca;
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = 'Turbulence Intensity [% [0, 1]]';
title('Turbulence Intensity');

hold off;

%% 7.2 Power Spectral Densities (Auto-Spectra)
%
% Other turbulence parameters include the TKE power- and cross-spectral 
% densities (i.e the power spectra), turbulent kinetic energy (TKE, i.e. 
% the variances of velocity vector components), Reynolds stress vector 
% (i.e. the co-variances of velocity vector components), TKE dissipation 
% rate, and TKE production rate. These quantities are primarily used to 
% inform and verify hydrodynamic and coastal models, which take some or 
% all of these quantities as input.
%
% The TKE production rate is the rate at which kinetic energy (KE) 
% transitions from a useful state (able to do "work" in the physics sense) 
% to turbulent; TKE is the actual amount of turbulent KE in the water; and 
% TKE dissipation rate is the rate at which turbulent KE is lost to 
% non-motion forms of energy (heat, sound, etc) due to viscosity. The power 
% spectra are used to depict and quantify this energy in the frequency 
% domain, and creating them are the first step in turbulence analysis.
%
% This section examines the power spectra, specifically the auto-spectra
% from the vertical beam. This can be done using the
% |calculate_power_spectral_density| function. The following code creates spectra at the
% middle water column, at a depth of 5 m, and uses a number of FFT's equal
% to 1/3 the bin size.

rng = 5;  % m - depth for analysis

[vel_up, vel_coords] = dolfyn_select(ds.vel_b5, 'range_b5', rng, 'method', 'nearest');

[U_data, u_coords] = dolfyn_select(ds_avg.U_mag, 'range', rng, 'method', 'nearest');

ds_avg = power_spectral_density(ds_avg, vel_up, ...
    'freq_units', 'Hz', ...
    'n_fft', floor(ds_avg.attrs.n_bin / 2), ...
    'field_name', 'auto_spectra_5m');

U = U_data;

% Create figure for PSD plot with velocity-colored spectra
fig = figure;
fig.Position = [100 100 500 500];  

% Set default font properties to match Python
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultAxesFontName', 'Times New Roman');

% Get PSD and frequency data
psd_data = ds_avg.auto_spectra_5m.data;  % [time_bins x freq_bins]
freq_data = ds_avg.auto_spectra_5m.coords.freq;

% Get velocity data for coloring (U was extracted earlier)
U_values = U_data(:);  % Make sure it's a column vector
U_max = max(U_values);

% Average spectra into 0.1 m/s velocity bins (match Python exactly)
% Python: np.arange(0.5, U_max, 0.1)
speed_bins = 0.5:0.1:U_max;
if speed_bins(end) ~= U_max
    speed_bins = [speed_bins, U_max];  % Ensure U_max is included
end
n_bins = length(speed_bins) - 1;

% Initialize storage for binned spectra
S_binned = zeros(n_bins, length(freq_data));
count_binned = zeros(n_bins, 1);

% Bin the spectra by velocity (match Python's groupby_bins behavior)
for i = 1:n_bins
    if i < n_bins
        % For all bins except last: [bin_start, bin_end)
        bin_mask = (U_values >= speed_bins(i)) & (U_values < speed_bins(i+1));
    else
        % For last bin: [bin_start, bin_end]
        bin_mask = (U_values >= speed_bins(i)) & (U_values <= speed_bins(i+1));
    end
    
    if sum(bin_mask) > 0
        S_binned(i, :) = mean(psd_data(bin_mask, :), 1, 'omitnan');
        count_binned(i) = sum(bin_mask);
    else
        S_binned(i, :) = NaN;
    end
end

% Create colormap (jet is similar Python's turbo)
cmap = jet(256);


% Create colors for each speed bin (match Python's BoundaryNorm)
cbar_max = 2.0;  % Match Python's cbar_max=2.0 exactly
bounds = 0.5:0.1:cbar_max;  % Match Python's np.arange(0.5, cbar_max, 0.1)

% Map speed_bins to colors using BoundaryNorm approach
speed_bin_centers = speed_bins(1:end-1) + diff(speed_bins)/2;
colors = zeros(n_bins, 3);
for i = 1:n_bins
    % Find which boundary interval this speed falls into
    bound_idx = find(bounds <= speed_bin_centers(i), 1, 'last');
    if isempty(bound_idx)
        bound_idx = 1;
    elseif bound_idx >= length(bounds)
        bound_idx = length(bounds) - 1;
    end
    
    % Map to colormap index
    cmap_idx = round((bound_idx - 1) / (length(bounds) - 1) * (size(cmap, 1) - 1)) + 1;
    cmap_idx = max(1, min(size(cmap, 1), cmap_idx));
    colors(i, :) = cmap(cmap_idx, :);
end

% Create main plot axes (match Python's subplots_adjust layout)
ax = axes('Position', [0.2 0.1 0.55 0.85]);  % left=0.2, right=0.75, bottom=0.1, top=0.95
hold on;

% Plot each velocity bin with its color
for i = 1:n_bins
    if count_binned(i) > 0 && ~all(isnan(S_binned(i, :)))
        loglog(freq_data, S_binned(i, :), 'Color', colors(i, :), 'LineWidth', 1.0);
    end
end

% Add -5/3 slope reference line
m = -5/3;
x = logspace(-1, 0.5, 50);  % Match Python's np.logspace(-1, 0.5)
y = 10^(-3) * x.^m;         % Match Python's 10**(-3) * x**m
loglog(x, y, '--', 'Color', 'black', 'LineWidth', 2);

% Add grid (match Python)
grid on;

% Specify logarithmic scaling
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');

% Set axis properties (match Python exactly)
set(ax, 'FontSize', 18, 'FontName', 'Times New Roman');
xlabel('Frequency [Hz]', 'FontSize', 18);
ylabel('PSD [m^2 s^{-2} Hz^{-1}]', 'FontSize', 18);  % Match Python's format
xlim([0.01, 1]);
ylim([0.0005, 0.1]);

% Add legend (match Python's formatting)
legend('f^{-5/3}', 'Location', 'northeast');

% Create colorbar (match Python's ColorbarBase implementation)
cax = axes('Position', [0.8 0.07 0.03 0.88]);  % Match Python's [0.8, 0.07, 0.03, 0.88]
% Create colorbar data
cbar_data = linspace(0.5, cbar_max, 256)';
imagesc(cax, [1], cbar_data, repmat(cbar_data, 1, 1));
colormap(cax, cmap);
% Format colorbar (match Python's formatting)
set(cax, 'YDir', 'normal', 'XTick', []);
set(cax, 'FontSize', 18, 'FontName', 'Times New Roman');
set(cax, 'YLim', [0.5, cbar_max]);
% Set colorbar ticks to match Python's bounds
yticks(cax, bounds);
yticklabels(cax, arrayfun(@(x) sprintf('%.1f', x), bounds, 'UniformOutput', false));

% Move y-axis (ticks and label) to the right side
set(cax, 'YAxisLocation', 'right');
ylabel(cax, 'Velocity [m/s]', 'FontSize', 18);

hold off;

%% 7.3 Instrument Noise
%
% The next step calculates the instrument's Doppler noise
% floor from the spectrum calculated above. (This analysis assumes that the noise floor of the vertical beam is the same as the
% noise floor of the other 4 beams). This provides a timeseries of the
% noise floor, which varies by instrument and with flow speed, at that
% depth bin.
%
% This section uses the |calculate_doppler_noise_level| function. The
% two inputs for this function are the power spectra and "pct_fN", the
% percent of the Nyquist frequency that the noise floor exists. Because in
% this particular dataset the noise floor is not visible, this section uses
% 90% or pct_fN=0.9 as an example. If the noise floor began at 0.4 Hz and
% ran to the maximum frequency of 0.5 Hz, the calculation would use pct_fN = 0.4 Hz / 0.5 Hz = 0.8.

ds_avg = calculate_doppler_noise_level(ds_avg, ...
    'psd_field', 'auto_spectra_5m', ...
    'pct_fN', 0.9, ...
    'field_name', 'noise_5m');

% Display noise statistics
noise_data = ds_avg.noise_5m.data;
doppler_noise_stats = table(...
    mean(noise_data(:), 'omitnan'), ...
    std(noise_data(:), 'omitnan'), ...
    min(noise_data(:)), ...
    max(noise_data(:)), ...
    'VariableNames', {'Mean_m_s', 'StdDev_m_s', 'Min_m_s', 'Max_m_s'});
doppler_noise_stats.Properties.Description = 'Doppler Noise Level Statistics';
disp(doppler_noise_stats);

%% 7.4 TKE Dissipation Rate
%
% Because the isotropic turbulence cascade (0.2 - 0.5 Hz) is visible at this
% depth bin (5 m altitude), the TKE dissipation rate can be calculated at this 
% location from the spectra itself. This can be done using 
% |calculate_dissipation_rate_LT83|, whose inputs are the power spectra, the 
% ensemble speed, the frequency range of the isotropic cascade, and the 
% instrument's noise.

% Frequency range of isotropic turbulence cascade
f_rng = [0.2, 0.5];  % Hz

% Create temporary 1D U_mag field for dissipation calculation  
% (since the PSD is 2D [time, freq], a 1D U_mag [time] is required)
ds_avg.U_mag_5m = struct();
ds_avg.U_mag_5m.data = U_data;  % Use the selected 1D data
ds_avg.U_mag_5m.dims = {'time'};  % 1D time dimension
ds_avg.U_mag_5m.coords = struct();  % Empty coords structure

% Calculate dissipation rate
ds_avg = calculate_dissipation_rate_LT83(ds_avg, ...
    'psd_field', 'auto_spectra_5m', ...
    'U_mag_field', 'U_mag_5m', ...
    'freq_range', f_rng, ...
    'noise', ds_avg.noise_5m, ...
    'field_name', 'dissipation_rate_5m');

% Display dissipation statistics
eps_data = ds_avg.dissipation_rate_5m.data;
dissipation_stats = table(...
    mean(eps_data(:), 'omitnan'), ...
    sum(~isnan(eps_data(:))), ...
    numel(eps_data), ...
    100 * sum(~isnan(eps_data(:))) / numel(eps_data), ...
    'VariableNames', {'Mean_m2_s3', 'Valid_Points', 'Total_Points', 'Valid_Percent'});
dissipation_stats.Properties.Description = 'TKE Dissipation Rate Statistics (5m depth)';
disp(dissipation_stats);

%% Calculate full profile of spectra and dissipation rates
%
% The previous section found the spectra and dissipation rate from a single depth
% bin at an altitude of 5 m from the seafloor, but typically analysis requires the
% spectra and dissipation rates from the entire measurement profile.

ds_avg = calculate_dissipation_rate_profile(ds_avg, ds, 'freq_range', f_rng);


%% Quality control for dissipation rate
%
% With a profile timeseries of dissipation rate calculated, the next step applies
% quality control (QC). Since individual spectra cannot be manually inspected
% to verify the isotropic turbulence cascade, QC becomes necessary for the output
% from |calculate_dissipation_rate_LT83| to verify that calculated values
% actually fall on a $f^{(-5/3)}$ slope. The function |check_turbulence_cascade_slope|
% provides this capability, using linear regression on the log-transformed LT83 equation
% (ref. to Lumley and Terray, 1983, see docstring) to calculate the spectral slope for the
% given frequency range.
%
% This analysis calculates the slope of each spectrum between 0.2 and
% 0.5 Hz. The following uses a cutoff of 20% for the error, but this can be lowered
% if erroneous estimations appear during visual inspection
% of the spectra.

% Display QC results from the profile function
valid_estimates = sum(ds_avg.qc_mask.data(:));
total_estimates = numel(ds_avg.qc_mask.data);
qc_results = table(...
    valid_estimates, ...
    total_estimates, ...
    100*valid_estimates/total_estimates, ...
    'VariableNames', {'Valid_Estimates', 'Total_Estimates', 'Pass_Rate_Percent'});
qc_results.Properties.Description = 'Quality Control Results for Dissipation Rate';
disp(qc_results);

% Apply quality control mask
threshold_for_difference_from_expected_slope = 0.20;  % 20%
slope_data = ds_avg.qc_slope.data;
mask = abs((slope_data - (-5/3)) / (-5/3)) <= threshold_for_difference_from_expected_slope;

% Apply mask to dissipation rate (match Python exactly)
ds_avg.dissipation_rate.data(~mask) = NaN;

%% Physical interpretation and ADCP limitations
%
% For this example dataset collected at 1 Hz, plotting the dissipation rate in a colormap shows
% significant missing data in the profile map. This 1 Hz sampling rate doesn't provide
% sufficient information for reliable dissipation rate estimations, and turbulence measurements
% approach the limits of ADCP measurement capabilities.
%
% Also, 1x10^-4 to 3x10^-4 $m^2/s^3$ is reasonable for a dissipation rate 
% estimate for the 1 - 1.5 m/s current speeds measured here. They can be a 
% magnitude greater for faster flow speeds, typically increase closer to the 
% seafloor, and depend heavily on bathymetry and regional hydrodynamics.

% Create dissipation rate visualization
fig = figure;
fig.Position = [100 100 800 400];

% Use dissipation rate profile from the function
dissipation_full = ds_avg.dissipation_rate.data;

time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
range = ds.coords.range;
depth = ds_avg.depth.data;

[T, R] = meshgrid(time, range);
pcolor(T, R, dissipation_full);
shading flat;
% colormap('turbo');
colormap(cmocean('thermal', 256));

hold on;
surface_elevation = plot(squeeze(time), squeeze(depth), 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'northeast');

xlabel('time b5');
ylabel('range');
ylim([0, max(depth) + 1]);

ax = gca;
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = 'TKE Dissipation Rate [m2 s-3]';
title('TKE Dissipation Rate');

hold off;

%% 7.5 Noise-Corrected Turbulence Intensity
%
% With the noise floor calculated for each ping, TI can be recalculated with instrument noise removed by using the |calculate_turbulence_intensity|
% function. Subtracting this from the non-noise corrected function
% shows a large difference at slower flow speeds, with an average
% difference of about 0.008 (0.8%). This approach also removes measurements
% where noise is high.


ds_avg = calculate_turbulence_intensity(ds_avg, ...
    'noise', ds_avg.noise_5m, ...
    'field_name', 'turbulence_intensity');

% Calculate and plot the difference
ti_diff = ds_avg.TI.data - ds_avg.turbulence_intensity.data;

% Create visualization of TI difference
fig = figure;
fig.Position = [100 100 800 400];

time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
range = ds_avg.coords.range;
depth = ds_avg.depth.data;

[T, R] = meshgrid(time, range);
pcolor(T, R, ti_diff');
shading flat;
colormap(cmocean('algae', 256));  % Green colormap similar to Python

hold on;

% Plot water surface depth
surface_elevation = plot(time, depth, 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'southeast');

xlabel('Time');
ylabel('range');
ylim([0, max(depth) + 1]);

ax = gca;
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = 'TI Difference';
title('TI Difference');

hold off;

%% 7.6 Reynolds Stress Components
%
% The next parameters to calculate are the Reynolds normal and shear 
% stresses (-$\overline{u_iu_j}$). Using the vertical beam on
% the ADCP allows calculation of the vertical TKE component from the
% along-beam velocity using the |calculate_turbulent_kinetic_energy|
% function. This function calculates TKE for any along-beam
% velocity.
%
% The so-called "beam-variance" equations can also estimate the 
% Reynolds stress tensor components (i.e. $\overline{u'}^2$, $\overline{v'}^2$, 
% $\overline{w'}^2$, $\overline{u'v'}$, $\overline{u'w'}^2$, $\overline{v'w'}^2$), 
% which define the normal and shear stresses acting on an element of water. 
% These equations are built into the functions |calculate_reynolds_stress_5beam| 
% and |calculate_reynolds_stress_4beam|.
%
% Both of these functions will give comparable results, but 
% |calculate_reynolds_stress_5beam| takes into account instrument tilt, and 
% |calculate_reynolds_stress_4beam| does not. Both will throw a tilt warning 
% if tilt is greater than 5 degrees.
%
% *5-beam ADCP considerations:*
%
% There are a couple caveats to calculating Reynolds stress tensor components:
%
% # Because this instrument only has 5 beams, only 5 of the 6 components can be found (6 unknowns, 5 knowns)
% # Because the ADCP's instrument (XYZ) axes weren't aligned with the flow during deployment, the direction these components are aligned to is unknown (i.e. the 'u' direction is not necessarily the streamwise direction)
% # It is possible to rotate the tensor, but all 6 components would be required to do so properly ("coupled ADCPs")
% # Measurements close to the seafloor can be suspect due to increased vertical flow. ADCPs operate under the "assumption of homogeneity", which means that they can only accurately measure consistent horizontal currents with relatively little vertical motion.
%
% The following section calculates the Reynolds stresses using the 
% |calculate_reynolds_stress_5beam| function, which calculates the individual 
% Reynolds stress tensor components and takes the same inputs: the raw dataset 
% in "beam" coordinates, the instrument Doppler noise, the ADCP's orientation 
% and its beam angle. It outputs both the normal stress ("tke_vec") and the 
% shear stress ("stress_vec") vector. Note, this function will drop at least 
% one warning every time it's run, primarily the coordinate system warning. 
% This function also requires the input raw data to be in beam coordinates, 
% so this section creates a copy of the raw data and rotates it to 'beam'. If not provided,
% this function will perform the rotation automatically.

ds_beam = rotate2(ds, 'beam');
% Extract mean noise value for 5-beam calculation
noise_mean = mean(ds_avg.noise_5m.data(:), 'omitnan');
ds_reynolds = calculate_reynolds_stress_5beam(ds_beam, ...
    'noise', noise_mean, ...
    'orientation', 'up', ...
    'beam_angle', 25, ...
    'align_with_shear', true);

% Merge Reynolds stress results into ds_avg (preserve existing fields like dissipation_rate)
ds_avg.tke_vec = ds_reynolds.tke_vec;
ds_avg.stress_vec = ds_reynolds.stress_vec;

% Display Reynolds stress statistics
tke_data = ds_avg.tke_vec.data;
stress_data = ds_avg.stress_vec.data;
reynolds_stats = table(...
    mean(tke_data(:), 'omitnan'), ...
    sum(~isnan(tke_data(:))), ...
    numel(tke_data), ...
    mean(abs(stress_data(:)), 'omitnan'), ...
    sum(~isnan(stress_data(:))), ...
    numel(stress_data), ...
    'VariableNames', {'TKE_Mean_m2_s2', 'TKE_Valid_Points', 'TKE_Total_Points', ...
                      'Stress_Mean_m2_s2', 'Stress_Valid_Points', 'Stress_Total_Points'});
reynolds_stats.Properties.Description = 'TKE and Reynolds Stress Statistics';
disp(reynolds_stats);

%%
%
% There is one other important thing to note on Reynolds stress measurements 
% by ADCPs: the minimum turbulence length scale that the ADCP is capable of 
% measuring increases with range from the instrument. This means the 
% instrument is only capable of measuring the stress transported by larger 
% and larger turbulent structures as the beams travel farther and farther 
% from the instrument head. One of the benefits of calculating w'w' from the 
% vertical beam is that it isn't limited by this beam spread issue, though 
% on its own it may not be particularly useful.

%% 7.7 TKE Production
%
% Though it can't be found from this deployment, the following explains how to estimate
% the TKE shear production rate. This calculation assumes that buoyancy
% production is negligible and the environment is a well-mixed tidal channel. MHKiT-DOLfYN does not include
% a specific function for production, but provides all necessary
% variables.
%
% It is possible to estimate production rates from either an ADV or an ADCP 
% aligned with the flow direction (so "X" would align with the principal flow 
% direction). The following estimation for shear production rates takes into 
% account both along- and cross-stream shear:
%
% $$P = -\overline{u'w'}\frac{du}{dz} - \overline{v'w'}\frac{dv}{dz}$$
%
% Note that the signs can be complex but are important here. The previous calculations found the
% Reynolds shear stresses $-\overline{u'w'}$ and $-\overline{v'w'}$ using
% the Reynolds stress equations. If ADV data is available, those estimations
% are preferred because the ADV's point measurement does not have the assumptions
% that ADCP measurements have.
%
% The velocity shear components can be found from the aptly named functions 
% |dudz| and |dvdz| in ADPBinner. These functions, which are useful alone in 
% their own right, approximate the shear in the velocity vector between respective 
% depth bins. There is always correlation between velocity measurements in adjacent 
% depth bins, based on ADCP operation principles, which is why "approximate" is 
% also used here for velocity shear.
%
% The velocity shear functions operate on the raw velocity vector in the principal 
% reference frame and need to be ensemble-averaged here. This can be done by 
% nesting the |d*dz| function within the ADPBinner's |mean| function. With the 
% ensemble shear known, all components can be combined to calculate a production
% estimation. If using ADV data, take the mean again of the "range" dimension.


% Calculate velocity shear components
ds_avg = calculate_velocity_shear(ds_avg, ...
    'vel_field', 'vel', ...
    'diff_style', 'centered');

% Extract Reynolds shear stress components 
% The updated calculate_reynolds_stress_5beam function now handles dimension alignment
% internally, so stress components should match velocity shear dimensions

% Access stress components - data is [range x time x tau]
upwp = squeeze(ds_avg.stress_vec.data(:, :, 2));  % u'w' component (2nd tau component)
vpwp = squeeze(ds_avg.stress_vec.data(:, :, 3));  % v'w' component (3rd tau component)

% Extract velocity shear
dudz = ds_avg.dudz.data;
dvdz = ds_avg.dvdz.data;
dwdz = ds_avg.dwdz.data;

% Extract wpwp from TKE vector (wpwp is the 3rd component: upup, vpvp, wpwp)
if isfield(ds_avg, 'tke_vec') && size(ds_avg.tke_vec.data, 3) >= 3
    wpwp = squeeze(ds_avg.tke_vec.data(:, :, 3));  % w'w' component (3rd tke component)
else
    error('wpwp component not found in tke_vec data');
end


% Calculate TKE production using main horizontal terms (standard oceanographic approach)
% P = -u'w'*du/dz - v'w'*dv/dz
% Note: w'w'*dw/dz term often excluded due to measurement limitations and smaller contribution

P = -upwp .* dudz - vpwp .* dvdz;


%% 7.8 TKE Balance
%
% Plotting the production rates and the ratio of production rates to
% dissipation rates provides understanding of the TKE balance. Production
% should typically be greater than 0, though negative values can indicate
% uncertainty. In a well mixed coastal environment,
% production and dissipation are expected to be approximately equal. These
% production estimates are possibly high because the stress components
% are not aligned with the flow (4x10^-3 m²/s³ is quite large), but if this
% were not the case, the conclusion would be that TKE is produced (kinetic energy
% is lost to turbulence) but not dissipated (turbulent energy is lost to
% entropy) at this location.

% Remove estimations below 0
P(P <= 0) = NaN;

% Plot production rate
fig = figure;
fig.Position = [100 100 viz_width viz_height];

time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
% Use the shear-aligned range coordinate (26 bins) instead of original range (28 bins)
if isfield(ds_avg, 'dudz') && isfield(ds_avg.dudz, 'coords') && isfield(ds_avg.dudz.coords, 'range')
    range = ds_avg.dudz.coords.range;  % 26 range bins from velocity shear
else
    % Fallback: create aligned range by trimming original range to match P dimensions
    original_range = ds_avg.coords.range;  % 28 bins
    n_range_p = size(P, 1);  % 26 bins
    range_start = 2;  % Skip first bin (centered difference)
    range_end = range_start + n_range_p - 1;  % Skip last bin
    range = original_range(range_start:range_end);
end
depth = ds_avg.depth.data;

% Create meshgrid
[T, R] = meshgrid(time, range);

% Ensure P has the right orientation for pcolor
% pcolor expects Z to match meshgrid dimensions: [length(range) x length(time)]
if size(P, 1) == length(range) && size(P, 2) == length(time)
    % P is already [range x time]
    pcolor(T, R, P);
elseif size(P, 1) == length(time) && size(P, 2) == length(range)
    % P is [time x range], transpose it
    pcolor(T, R, P');
end

shading flat;
% colormap('turbo');
colormap(cmocean('thermal', 256));

hold on;
surface_elevation = plot(squeeze(time), squeeze(depth), 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'northeast');

xlabel('Time');
ylabel('Profile Range [m]');
ylim([0, max(depth) + 1]);

ax = gca;
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = 'stress\_vec';
title('TKE Production');

hold off;

% Plot difference between production and dissipation rate
fig = figure;
fig.Position = [100 100 viz_width viz_height];


% Extract dissipation rate from full profile (match Python: ds_avg["dissipation_rate"].values)
% Align dimensions: P is [26 x time], dissipation_rate is [28 x time]
% Use same range subset as production (centered difference removes first and last bins)
range_start = 2;  % Skip first bin (centered difference)
range_end = range_start + size(P, 1) - 1;  % Take 26 bins to match P
eps_profile = ds_avg.dissipation_rate.data(range_start:range_end, :);

% Calculate balance (match Python: P - ds_avg["dissipation_rate"].values)
balance = P - eps_profile;

pcolor(T, R, balance);
shading flat;
% Blue to Red
colormap(cmocean('balance', 256));

hold on;
surface_elevation = plot(squeeze(time), squeeze(depth), 'Color', surface_elevation_blue, 'LineWidth', 2);
legend(surface_elevation, 'Surface Elevation [m]', 'Location', 'northeast');

xlabel('Time');
ylabel("Profile Range [m]");
ylim([0, max(depth) + 1]);

ax = gca;
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = 'stress\_vec';
title('TKE Balance');

hold off;
