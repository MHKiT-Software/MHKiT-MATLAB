%% Analyzing ADCP Data with MHKiT
% The following example illustrates a straightforward workflow for analyzing 
% Acoustic Doppler Current Profiler (ADCP) data utilizing MHKiT. MHKiT has integrated 
% the Doppler Oceanographic Library for pYthoN (DOLfYN) codebase as a module to 
% facilitate ADCP and Acoustic Doppler Velocimetry (ADV) data processing.
% 
% Here is a standard workflow for ADCP data analysis:
% 
% *1. Import Data*
% 
% *2. Review, QC, and Prepare the Raw Data:*
%% 
% # Calculate or verify the correctness of depth bin locations
% # Discard data recorded above the water surface or below the seafloor
% # Assess the quality of velocity, beam amplitude, and/or beam correlation 
% data
% # Rotate Data Coordinate System
%% 
% *3. Data Averaging:*
%% 
% * If not already executed within the instrument, average the data into time 
% bins of a predetermined duration, typically between 5 and 10 minutes
%% 
% *4. Calculation of Speed and Direction*
% 
% *5. Visualizations*
% 
% *6. Saving and Loading DOLfYN datasets*
% 
% *7. Turbulence Statistics:*
%% 
% * Turbulence Intensity (TI)
% * Power Spectral Densities
% * Instrument Noise
% * Turbulent Kinetic Energy (TKE) Dissipation Rate
% * Noise-corrected TI
% * TKE Componenets
% * TKE Production
% * TKE Balance
%% 1. Importing Raw Instrument Data
% One of DOLfYN's key features is its ability to directly import raw data from 
% an Acoustic Doppler Current Profiler (ADCP) right after it has been transferred. 
% In this instance, we are using a Nortek Signature1000 ADCP, with the data stored 
% in files with an '.ad2cp' extension. This specific dataset represents several 
% hours of velocity data, captured at 1 Hz by an ADCP mounted on a bottom lander 
% within a tidal inlet. The list of instruments compatible with DOLfYN can be 
% found in the <https://mhkit-software.github.io/MHKiT/mhkit-python/api.dolfyn.html 
% MHKiT DOLfYN documentation>.
% 
% We'll start by importing the raw data file downloaded from the instrument. 
% The |dolfyn_read| function processes the input raw file and converts the data 
% into a MATLAB struct, typically called a Dataset, or |ds| for short. This Dataset 
% includes several groups of variables:
%% 
% # *Velocity*: Recorded in the coordinate system saved by the instrument (beam, 
% XYZ, ENU)
% # *Beam Data*: Includes amplitude and correlation data
% # *Instrumental & Environmental Measurements*: Captures the instrument's bearing 
% and environmental conditions
% # *Orientation Matrices*: Used by DOLfYN for rotating through different coordinate 
% frames.
%% 
% The |dolfyn_read| function parses the raw ACDP file into MATLAB.

tic;
ds = dolfyn_read('./examples/data/dolfyn/Sig1000_tidal.ad2cp');
toc
%% 
% To see what is in the data set explore the |ds| variable either in the workspace 
% or below

ds
%% 
% The function |dolfyn_plot| provides a simple api to visualize the data read 
% by DOLFyN. Pass in the dataset valiable, the dimension name and index. |dolfyn_plot| 
% uses data found within the dataset attributes to determine the plot type and 
% labels. In the following code section we visualize the velocity of the East, 
% North, and U1 directions of the ADCP data over time and depth.

figure('Position', [100 100 1200 300]);

% Create three subplots in a 1x3 grid
subplot(1,3,1)
dolfyn_plot(ds,'vel','dim',1)

subplot(1,3,2)
dolfyn_plot(ds,'vel','dim',2)

subplot(1,3,3)
dolfyn_plot(ds,'vel','dim',3)

% Adjust the spacing between subplots
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1 0.1 0.8 0.3])
%% 2. Initial Steps for Data Quality Control (QC)
% *2.1: Set the Deployment Height*
% 
% When using Nortek instruments, the deployment software does not factor in 
% the deployment height. The deployment height represents the position of the 
% Acoustic Doppler Current Profiler (ADCP) within the water column.
% 
% In this context, the center of the first depth bin is situated at a distance 
% that is the sum of three elements:
%% 
% # Deployment height (the ADCP's position in the water column)
% # Blanking distance (the minimum distance from the ADCP to the first measurement 
% point)
% # Cell size (the vertical distance of each measurement bin in the water column)
%% 
% To ensure accurate readings, it is critical to calibrate the 'range' coordinate 
% to make '0' correspond to the seafloor. This calibration can be achieved using 
% the |set_range_offset| function. This function is also useful when working with 
% a down-facing instrument as it helps account for the depth below the water surface.
% 
% For those using a Teledyne RDI ADCP, the TRDI deployment software will prompt 
% the user to specify the deployment height/depth during setup. If there's a need 
% for calibration post-deployment, the |set_range_offset| function can be utilized 
% in the same way as described above.
% 
% The ADCP transducers were measured to be 0.6 meters from the feet of the lander

ds = set_range_offset(ds, 0.6);
%% 
% We ran verify the range values were updated correctly by printing the min, 
% max and spacing (diff) of the range values defined in |ds.coords.range|

fprintf('Range Max: %.2f\n', max(ds.coords.range))
fprintf('Range Min: %.2f\n', min(ds.coords.range))
fprintf('Range Spacing Mean: %.2f\n', mean(diff(ds.coords.range)))
fprintf('Range Spacing Max: %.2f\n', max(diff(ds.coords.range)))
fprintf('Range Spacing Min: %.2f\n', min(diff(ds.coords.range)))
%% 
% *2.2. Discard Data Above Surface Level*
% 
% To reduce computational load, we can exclude all data at or above the water 
% surface level. Since the instrument was oriented upwards, we can utilize the 
% pressure sensor data along with the function |water_depth_from_pressure|. However, 
% this approach necessitates that the pressure sensor was calibrated or 'zeroed' 
% prior to deployment. If the instrument is facing downwards or doesn't include 
% pressure data, the function |water_depth_from_pressure| can be used to detect 
% the seabed or water surface.
% 
% It's important to note that Acoustic Doppler Current Profilers (ADCPs) do 
% not measure water salinity, so the user will need to supply this information 
% to the function. The dataset returned by this function includes an additional 
% variable, "depth". If |water_depth_from_pressure| is invoked after |set_range_offset|, 
% "depth" represents the distance from the water surface to the seafloor. Otherwise, 
% it indicates the distance to the ADCP pressure sensor. The |water_depth_from_pressure| 
% functions allows the user to specified salinity value in Practical Salinity 
% Units (psu).

water_salinity_psu = 31;
ds = water_depth_from_pressure(ds, 'salinity', water_salinity_psu);
%% 
% Visualizing Depth

figure;
dolfyn_plot(ds, 'depth');
title("Distance from Water Surface to Seafloor [m]");
%% 
% After determining the "depth", the |remove_surface_interference| function 
% can be used to discard data in depth bins near or above the actual water surface. 
% This function calculates a range limit based on the beam angle and cell size, 
% where surface interference is expected to occur at distances greater than |range 
% * cos(beam angle) - cell size|. The beam angle accounts for the acoustic spread 
% of the ADCP signal (typically 20-25 degrees). To exclude the found surface interferance 
% data DOLFyN set and sets data above this threshold to NaN.

ds = remove_surface_interference(ds);
%% 
% Visualizing removal of surface interference

figure('Position', [100 100 1200 300])

subplot(1,3,1)
dolfyn_plot(ds,'vel','dim',1)

subplot(1,3,2)
dolfyn_plot(ds,'vel','dim',2)

subplot(1,3,3)
dolfyn_plot(ds,'vel','dim',3)

% Adjust the spacing between subplots
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1 0.1 0.8 0.3])
%% 
% *Correlation*
% 
% It's beneficial to also review data from the other beams. A significant portion 
% of this data is of high quality. To avoid discarding valuable data with lower 
% correlations, which could be due to natural variations, we can use the |correlation_filter|. 
% This function assigns a value of NaN (not a number) to velocity values corresponding 
% to correlations below 50%.
% 
% However, it's important to note that the correlation threshold is dependent 
% on the specifics of the deployment environment and the instrument used. It's 
% not unusual to set a threshold as low as 30%, or even to forgo the use of this 
% function entirely.

ds = correlation_filter(ds, "thresh", 50);
%% 
% Visualizing Correlation Data

figure('Position', [100 100 1200 300])

subplot(1,4,1)
dolfyn_plot(ds,'corr','dim',1)

subplot(1,4,2)
dolfyn_plot(ds,'corr','dim',2)

subplot(1,4,3)
dolfyn_plot(ds,'corr','dim',3)

subplot(1,4,4)
dolfyn_plot(ds,'corr','dim',4)

% Adjust the spacing between subplots
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1 0.1 0.8 0.3])
%% 
% *2.4 Rotate Data Coordinate System*
% 
% After cleaning the data, the next step is to rotate the velocity data into 
% accurate East, North, Up (ENU) coordinates.
% 
% ADCPs utilize an internal compass or magnetometer to determine magnetic ENU 
% directions. You can use the set_declination function to adjust the velocity 
% data according to the magnetic declination specific to your geographical coordinates. 
% This declination can be looked up online for specific coordinates.
% 
% Instruments save vector data in the coordinate system defined in the deployment 
% configuration file. To make this data meaningful, it must be transformed through 
% various coordinate systems ("beam"<->"inst"<->"earth"<->"principal"). This transformation 
% is accomplished using the |rotate2| function. If the "earth" (ENU) coordinate 
% system is specified, DOLfYN will automatically rotate the dataset through the 
% required coordinate systems to reach the "earth" coordinates.
% 
% In this case, since the ADCP data is already in the "earth" coordinate system, 
% the |rotate2| function will return the input dataset without modifications. 
% The |set_declination| function will work no matter the coordinate system.

ds = set_declination(ds, 15.8);
ds = rotate2(ds, 'earth');
%% 
% To rotate into the principal frame of reference (streamwise, cross-stream, 
% vertical), if desired, we must first calculate the depth-averaged principal 
% flow heading and add it to the dataset attributes. Then the dataset can be rotated 
% using the same |rotate2| function.

ds.attrs.principal_heading = calc_principal_heading(squeeze(ds.vel.data), true);
disp(ds.attrs.principal_heading);
ds_streamwise = rotate2(ds, 'principal');
%% 
% Visualize Streamwise Velocity

figure('Position', [100 100 1200 300])

% Create three subplots in a 1x3 grid
subplot(1,3,1)
dolfyn_plot(ds_streamwise,'vel','dim',1)

subplot(1,3,2)
dolfyn_plot(ds_streamwise,'vel','dim',2)

subplot(1,3,3)
dolfyn_plot(ds_streamwise,'vel','dim',3)

% Adjust the spacing between subplots
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1 0.1 0.8 0.3])
%% 3. Average the Data
% As this deployment was configured in "burst mode", a standard step in the 
% analysis process is to average the velocity data into time bins.
% 
% However, if the instrument was set up in an "averaging mode" (where a specific 
% profile and/or average interval was set, for instance, averaging 5 minutes of 
% data every 30 minutes), this step would have been performed within the ADCP 
% during deployment and can thus be skipped.
% 
% The function |average_by_dimension| is used to average the data by dimension. 
% To average the data into time bins (also known as ensembles), the sampling rate 
% and averaging duration must be defined. DOLfYN stores the sampling frequency 
% [Hz] in |.attrs.fs|. The averaging duration should be specified in seconds. 
% We can then use these values to calculate the number of samples to average.

sampling_frequency_hz = ds.attrs.fs;
averaging_frequency_seconds = 300; % 5 minutes
averaging_samples = int32(averaging_frequency_seconds / sampling_frequency_hz);
%% 
% Once the averaging samples count has been calculated we can average the dataset 
% by time by passing the original dataset, the number of samples to average, and 
% the dimension (|time|).

ds_avg = average_by_dimension(ds, averaging_samples, 'time');
%% 
% Plotting the Summary Data

figure('Position', [100 100 1200 300])  % Create a wide figure to accommodate the row

% Create three subplots in a 1x3 grid
subplot(1,3,1)
dolfyn_plot(ds_avg,'vel','dim',1)

subplot(1,3,2)
dolfyn_plot(ds_avg,'vel','dim',2)

subplot(1,3,3)
dolfyn_plot(ds_avg,'vel','dim',3)

% Adjust the spacing between subplots
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1 0.1 0.8 0.3])
%% 4. Calculating Current Speed and Direction
% The |calculate_horizontal_velocity_and_direction| function uses the ADCP east 
% and north vectors to compute current and direction. The variables |U_mag| and 
% |U_dir| will be added to the output dataset. |U_dir| is the "to" direction conversion 
% and direction is in degrees clockwise from true north. Note: ACDP flow direction 
% uses the "to" direction convention. Flow direction is the direction the current 
% is flowing towards. Wind direction convention is the direction the current is 
% coming from, and waves (typically wind driven) use the same "from" direction 
% convention.

ds_avg = calculate_horizontal_velocity_and_direction(ds_avg);

disp(ds_avg.U_mag);
disp(ds_avg.U_dir);
%% 5. Visualizing Current Speed and Direction
% Visualize Current Speed

fig = figure;
fig.Position = [100 100 800 500];

ax = axes('Position', [0.14 0.14 0.8 0.74]);

time = datetime(ds_avg.coords.time, 'ConvertFrom', 'posixtime');
range = ds_avg.coords.range;
speed = squeeze(ds_avg.U_mag.data);
speed = speed';
depth = ds_avg.depth.data;

[T, R] = meshgrid(time, range);
pcolor(T, R, speed);
shading flat;

blues = blues_colormap(256);
colormap(blues);

set(gca, 'Layer', 'top', 'Box', 'on')
% Force redraw of axes
set(gca, 'XGrid', 'off', 'YGrid', 'off')
set(gca, 'TickDir', 'out')
set(gca, 'LineWidth', 1)

hold on;

% Add a line denoting water surface depth using the darkest blue from the colormap
darkest_blue = blues(end,:);
plot(time, depth, 'Color', darkest_blue, 'LineWidth', 2);

% Set axis labels and limits
xlabel('Time');
ylabel('Altitude [m]');
ylim([0 12]);

% Format time
ax.XAxis.TickLabelFormat = 'HH:mm';

c = colorbar;
c.Label.String = ds_avg.U_mag.long_name + " [" + ds_avg.U_mag.units + "]";
title("Current Speed [m] and Surface Elevation [m]")

% Turn hold off
hold off;
%% 
% Visualize Current "To" Direction (Degrees CW from True North)

fig = figure;
fig.Position = [100 100 800 500];

ax = axes('Position', [0.14 0.14 0.8 0.74]);

direction = squeeze(ds_avg.U_dir.data);
direction = direction';

pcolor(T, R, direction);
shading flat;
colormap(twilight);

set(gca, 'Layer', 'top', 'Box', 'on')
% Force redraw of axes
set(gca, 'XGrid', 'off', 'YGrid', 'off')
set(gca, 'TickDir', 'out')
set(gca, 'LineWidth', 1)

hold on;

% Plot water surface depth
plot(time, depth, 'Color', darkest_blue, 'LineWidth', 2);

xlabel('Time');
ylabel('Altitude [m]');
ylim([0 12]);

ax.XAxis.TickLabelFormat = 'HH:mm';


c = colorbar;
c.Label.String = ds_avg.U_dir.long_name + " [deg CW from true N]";

title("Horizontal Velocity 'To' Direction [deg] and Surface Elevation [m]")

% Turn hold off
hold off;
%% 6. Saving Data
% Datasets can be saved and loaded using the write_netcdf and read_netcdf functions, 
% respectively.

% Uncomment these lines to save and load to your current working directory
% write_netcdf(ds, 'your_data.nc');
% ds_saved = read_netcdf('your_data.nc');
%% Upcoming Features
% In the next release of MHKiT-MATLAB turbulence statistics calculation and 
% visualization will be added