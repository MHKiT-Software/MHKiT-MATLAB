%% MHKiT WEC-Sim Power Performance Example
% This notebook demonstrates using the MHKiT wave module and WEC-Sim together 
% to perform a resource characterization study in MHKiT, simulate representative 
% cases with WEC-Sim, and visualize the results in MHKiT to estimate MAEP (Mean 
% Annual Energy Production).
% 
% 1. Characterize the available resource at a location
% 
% - See the PacWave example notebook
% 
% 2. Write a WEC-Sim batch file for the given clusters
% 
% 3. Simulate the device _in WEC-Sim_.
% 
% - Ensure that the spectra used in WEC-Sim is identical to the one used in 
% MHKiT.
% 
% 4. Load WEC-Sim batch results
% 
% 5. Assess results and visualize quantities of interest
% 
% This example uses WEC-Sim to simulate the [Oscillating Surge Wave Energy Converter 
% (OSWEC)](https://wec-sim.github.io/WEC-Sim/main/user/tutorials.html#oscillating-surge-wec-oswec), 
% a flap-type device.
%% 1. Characterize the available resource at a location
% This example will use an abbreviated version of `PacWave_resource_characterization_example.ipynb`. 
% 
% For full details on downloading, calculating, and visualizing the k-means 
% clusters representation of the site's wave resouce, see that example.
% 
% We select the N=32 cluster as it's total energy flux is closet to the total 
% energy flux of the site considering all wave conditions. We will load the PacWave 
% example output, which can be easily saved after running the example with the 
% command `results[32].to_csv("pacwave_cluster_32.csv")`. We will start this example 
% by reading in that csv output and formatting it for WEC-Sim.

% Relative location and filename of simulated WEC-Sim data (run with mooring)
filename = './data/wave/pacwave_cluster_32.csv'
results = readtable(filename)
%% 2. Write a WEC-Sim batch file for the given clusters
% WEC-Sim MCR (multiple condition run) files should contain a structure `mcr` 
% that contains two variables: `header` and `cases`. Each column of `header` and 
% `cases` denotes a variable and it's value respectively. Each row is distinct 
% simulation. WEC-Sim defines waves using the significant wave height and peak 
% period. We will isolate these values from the results of the cluster analysis 
% and create a dictionary that is written to the `.mat` file.

mcr = struct();
mcr.header = {'waves.height','waves.period'};
mcr.cases = [results.Hm0 results.Tp]
save('mcr_mhkit.mat', 'mcr');
%% 3. Simulate the device _in WEC-Sim_
% Now that the MCR file is created, we need to go simulate WEC performance in 
% these wave conditions using WEC-Sim. To recreate the data used in the next step, 
% use the created MCR file with WEC-Sim's OSWEC example. For an accurate comparison 
% to the power calculated in the resource characterization, we should ensure that 
% the WEC-Sim cases use irregular JONSWAP wave spectra as in the PacWave example.
% 
% 
% 
% For convenience in this demonstration, we enforce OSWEC model stability in 
% the extreme wave conditions by arbitrarily applying a large PTO stiffness and 
% damping in the wecSimInputFile.m:

% pto(1).stiffness = 1e5;
% pto(1).damping = 5e7;
%% 
% To reduce the amount of extranenous data saved for this example, we limit 
% the WEC-Sim output to the PTO's power output in the `userDefinedFunctions.m` 
% script:

% if exist('imcr','var')
%     if imcr == 1
%         nmcr = size(mcr.cases,1);
%         power = nan(1, nmcr);
%     end
% 
%     iRampEnd = simu.rampTime./simu.dtOut + 1;
%     power(imcr) = -mean(output.ptos(1).powerInternalMechanics(iRampEnd:end,5));
% 
%     if imcr == nmcr
%         % % Save output
%         save('mcr_mhkit_power.mat', 'power');
%     end
% end
% bodies = output.bodies;
%% 4. Load WEC-Sim batch results
% Note that in this example we do not save the entire WEC-Sim `output` structure 
% for each case. See the `wecsim_example.ipynb` for information on loading the 
% WEC-Sim responseClass. Here, the output is one array of average power output 
% that we will load and compare to the resource characterization.
% 
% Note that the power output [W] is significantly larger than the energy flux 
% [W/m] due to the large width of the OSWEC. 

% Relative location and filename of simulated WEC-Sim data (run with mooring)
filename = './data/wave/mcr_mhkit_power.mat';
% Load the WEC-Sim output data which contains the variable `power`.
load(filename)
results.Power = power'
%% 5. Assess results and visualize quantities of interest
% Now that we have loaded the OSWEC's modeled power, we can assess it's performance 
% relative to the incoming wave and calculate the mean annual energy production 
% (MAEP) using MHKiT.

results.CW = capture_length(results.Power, results.J)';

oswec_width = 18;
results.CWR = results.CW / oswec_width
%% 
% Calculate and display the mean annual energy production.

CW.values = results.CW;
J.values = results.J;
weights.values = results.weights;
MAEP = mean_annual_energy_production_matrix(CW, J, weights) / 1000 % kWh