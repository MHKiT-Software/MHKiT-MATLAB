function [clmat,maep_matrix] = power_performance_workflow(S, h, P, statistic, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     High-level function to compute power performance quantities of
%     interest following IEC TS 62600-100 Ed. 2.0 en 2024 for given wave
%     spectra.
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
%        Capture width statistics for plotting
%        options include: "mean", "std", "median",
%        "count", "sum", "min", "max", and "frequency".
%        Note that "std" uses a degree of freedom of N in accordance with
%        Formula D.5 of IEC TS 62600-100 Ed. 2.0 en 2024.
%        To output capture width matrices for multiple binning parameters,
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
%       Capture width matrix
%
%   maep_matrix: float
%       Mean annual energy production
%
% Examples
% --------
%     CDIP Example: Station Number 225: Kaneohe Bay, WETS, Oahu, HI
%     https://cdip.ucsd.edu/themes/cdip?pb=1&u2=s:225:st:1&d2=p70
%     >> station_number = '225';
%     >> data_type = 'realtime';
%     >> years = 2025;
%     >> parameters = {'waveEnergyDensity'};
%     >> data = cdip_request_parse_workflow('station_number', station_number, 'data_type', data_type, 'years', years, 'parameters', parameters);
%     >> % cast to double: fzero (inside wave_number, called by energy_flux
%     >> % below) requires double, but CDIP returns frequency as single.
%     >> S.frequency = double(data.metadata.wave.waveFrequency);  % [Hz], 64x1 double
%     >> % waveEnergyDensity is stored [time x frequency]; transpose to
%     >> % match MHKiT's frequency-as-rows, spectra-as-columns convention.
%     >> S.spectrum = data.data.wave2D.waveEnergyDensity';  % [m^2/Hz], 64x17520 double
%     >> S.time = data.data.wave.waveTime;  % 17520x1 datetime
%
%     >> % Water depth at CDIP 225 (WETS Berth 3), per PacIOOS buoy
%     >> % metadata: "Moored in water 80 meters deep".
%     >> h = 80;  % [m]
%
%     >> % SYNTHETIC power output - NOT measured. CDIP 225 is a wave buoy
%     >> % only, with no WEC device attached, so no real matching power
%     >> % record exists for this site. A real capture width matrix needs
%     >> % power telemetry from an actual device.
%     >> rng(1);
%     >> P = 50000 + 20000*randn(size(S.time));  % [W], illustrative only, 17520x1 double
%
%     >> [clmat, maep_matrix] = power_performance_workflow(S, h, P, "mean");
%     maep_matrix =
%         456328960.94  % [W*h]
%
%     >> clmat.mean
%     ans = 
%       struct with fields:
%           values: [9×14 double]
%             stat: 'mean'
%         Hm0_bins: [-0.2500 0.2500 0.7500 1.2500 1.7500 2.2500 2.7500 3.2500 3.7500]
%          Te_bins: [0.5000 1.5000 2.5000 3.5000 4.5000 5.5000 … ] (1×14 double)
%
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
    error('MHKiT:power_performance_workflow:InvalidInput', ...
        'S must be a structure with spectrum and frequency fields.');
end

if any([~isnumeric(h), ~isnumeric(P)])
    error('MHKiT:power_performance_workflow:InvalidInput', 'h and P must be numeric.');
end

% energy_period/significant_wave_height/energy_flux accept S.spectrum as
% a full NxK matrix directly, one column per spectrum - no per-spectrum
% loop needed.
if isnumeric(options.frequency_bins)
    Te = energy_period(S, options.frequency_bins);
    Hm0 = significant_wave_height(S, options.frequency_bins);
else
    Te = energy_period(S);
    Hm0 = significant_wave_height(S);
end
J = energy_flux(S, h, 'rho', options.rho, 'g', options.g);

% calculating capture width with power and wave flux in vectors
% Ensure P is a column vector to match J
P = P(:);
CW = capture_width(P, J);

% Need to set our Hm0 and Te bins for the capture width matrix
Hm0_bins = -0.5:0.5:max(fix(Hm0))+0.5; % Input is min, max, and n indecies for vector
Hm0_bins = Hm0_bins+0.25 ;
Te_bins = 0:1:max(fix(Te));
Te_bins = Te_bins+0.5;

% Calculate the necessary capture width matrices for each statistic based
% on IEC/TS 62600-100
clmat.mean   = capture_width_matrix(Hm0, Te, CW, "mean", Hm0_bins, Te_bins);
clmat.std    = capture_width_matrix(Hm0, Te, CW, "std", Hm0_bins, Te_bins);
clmat.median = capture_width_matrix(Hm0 ,Te, CW, "median", Hm0_bins, Te_bins);
clmat.count  = capture_width_matrix(Hm0 ,Te, CW, "count", Hm0_bins, Te_bins);
clmat.sum    = capture_width_matrix(Hm0 ,Te, CW, "sum", Hm0_bins, Te_bins);
clmat.min    = capture_width_matrix(Hm0 ,Te, CW, "min", Hm0_bins, Te_bins);
clmat.max    = capture_width_matrix(Hm0 ,Te, CW, "max", Hm0_bins, Te_bins);
clmat.freq   = capture_width_matrix(Hm0 ,Te, CW, "frequency", Hm0_bins, Te_bins);

% Create wave energy flux matrix using statistic
jmat = wave_energy_flux_matrix(Hm0, Te, J, "mean", Hm0_bins, Te_bins);

% Calcaulte MAEP from matrix
maep_matrix = mean_annual_energy_production_matrix(clmat.mean, jmat, clmat.freq);
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
        error('MHKiT:power_performance_workflow:InvalidInput', ...
            ['statistic must be a string or string array defined by one or ' ...
             'multiple of the following: "mean", "std", "median", "count", ' ...
             '"sum", "min", "max", "frequency".']);
    end
end

end