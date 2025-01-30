function [clmat,maep_matrix] = power_performance_workflow(S, h, P, statistic, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     High-level function to compute power performance quantities of
%     interest following IEC TS 62600-100 for given wave spectra.
%
% Parameters
% ------------
%   S: structure with fields:
%           S.spectrum: Spectral Density [m^2/Hz]
%           S.frequency: frequency [Hz]
%           S.time : time [datetime]
%
%   h: integer
%        Water depth [m]
%
%   P: array or vector
%        Power [W]
%
%   statistic: string or array of strings
%        Capture length statistics for plotting
%        options include: "mean", "std", "median",
%        "count", "sum", "min", "max", and "frequency".
%        Note that "std"
%        uses a degree of freedom of 1 in accordance with IEC/TS 62600-100.
%        To output capture length matrices for multiple binning parameters,
%        define as a string array: statistic = ["", "", ""];
%
%   savepath: string (optional)
%        Path to save figure.
%        to call: power_performance_wave(S,h,P,statistic,"savepath",savepath)
%
%   rho: float (optional)
%        Water density [kg/m^3]
%        to call: power_performance_wave(S,h,P,statistic,"rho",rho)
%
%   g: float (optional)
%        Gravitational acceleration [m/s^2]
%        to call: power_performance_wave(S,h,P,statistic,"g",g)
%
%   frequency_bins: vector (optional)
%      Bin widths for frequency of S. Required for unevenly sized bins
%
% Returns
% ---------
%   cl_matrix: figure
%       Capture length matrix
%
%   maep_matrix: float
%       Mean annual energy production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S
    h
    P
    statistic
    options.rho = 1025;
    options.g = 9.80665;
    options.frequency_bins = "";
    options.savepath = "";
end

if ~isstruct(S)
    ME = MException('MATLAB:power_performance_wave','S must be a structure. See API header for formatting');
    throw(ME);
end

if any([~isnumeric(h), ~isnumeric(P)])
    ME = MException('MATLAB:power_performance_wave','h and P must be numbers');
    throw(ME);
end

% Check dimensions of the spectrum data
[num_frequencies, num_spectra] = size(S.spectrum);

% Initialize arrays to store results
Te = zeros(num_spectra, 1);
Hm0 = zeros(num_spectra, 1);
J = zeros(num_spectra, 1);

% Process each spectrum separately
for i = 1:num_spectra
    % Create a temporary structure for the current spectrum
    S_temp = struct();
    S_temp.spectrum = S.spectrum(:,i);
    S_temp.frequency = S.frequency;

    % Compute the energy periods from spectra data
    Te(i) = energy_period(S_temp);

    % Compute the significant wave height from spectra data
    Hm0(i) = significant_wave_height(S_temp);

    % Compute the energy flux from spectra data and water depth
    J(i) = energy_flux(S_temp, h);
end

% calculating capture length with power and wave flux in vectors
L = capture_length(P,J);

% Need to set our Hm0 and Te bins for the capture length matrix
Hm0_bins = -0.5:0.5:max(fix(Hm0))+0.5; % Input is min, max, and n indecies for vector
Hm0_bins = Hm0_bins+0.25 ;
Te_bins = 0:1:max(fix(Te));
Te_bins = Te_bins+0.5;

% Calculate the necessary capture length matrices for each statistic based
% on IEC/TS 62600-100
clmat.mean = capture_length_matrix(Hm0,Te,L,"mean",Hm0_bins,Te_bins);
clmat.std = capture_length_matrix(Hm0,Te,L,"std",Hm0_bins,Te_bins);
clmat.median = capture_length_matrix(Hm0,Te,L,"median",Hm0_bins,Te_bins);
clmat.count = capture_length_matrix(Hm0,Te,L,"count",Hm0_bins,Te_bins);
clmat.sum = capture_length_matrix(Hm0,Te,L,"sum",Hm0_bins,Te_bins);
clmat.min = capture_length_matrix(Hm0,Te,L,"min",Hm0_bins,Te_bins);
clmat.max = capture_length_matrix(Hm0,Te,L,"max",Hm0_bins,Te_bins);
clmat.freq = capture_length_matrix(Hm0,Te,L,"frequency",Hm0_bins,Te_bins);

% Create wave energy flux matrix using statistic
jmat = wave_energy_flux_matrix(Hm0,Te,J,"mean",Hm0_bins,Te_bins);
% Calcaulte MAEP from matrix
maep_matrix = mean_annual_energy_production_matrix(clmat.mean,jmat,clmat.freq);
stats_cell = {'mean', 'std', 'median','count', 'sum', 'min', 'max','frequency'};

% Capture Length Matrix using statistic
cl_matrix = [];
len = strlength(options.savepath);
for i = 1:length(statistic)
    if any(strcmp(stats_cell,statistic(i)))
        figure('Name',sprintf('Capture Length Matrix %s', statistic(i)),'NumberTitle','off')
        cl_matrix(i) = plot_matrix(clmat.(statistic(i)),"Capture Length");
        name = [options.savepath, filesep, sprintf('Capture Length Matrix %s', statistic(i)), '.png'];

        if len > 1
            saveas(cl_matrix(i), name);
        end
    else
         ME = MException('MATLAB:power_performance_wave',...
             'statistic must be a string or string array defined', ...
             'by one or multiple of the following: "mean", "std", "median","count", "sum", "min", "max", "frequency"');
         throw(ME);
    end
end

end
