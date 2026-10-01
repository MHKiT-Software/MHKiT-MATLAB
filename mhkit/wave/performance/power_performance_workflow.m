function [cwmat, maep_matrix] = power_performance_workflow(S, h, P, statistic, options)

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
%   statistic: string or string array
%        Capture width statistics to compute and plot, one figure each.
%        Options: "mean", "std", "median", "count", "sum", "min", "max",
%        "probability", "frequency". "probability" and "frequency" are the
%        same statistic. "std" is the population standard deviation,
%        matching MHKiT-Python. "mean" and "probability" are always
%        computed because MAEP needs them.
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
%   cwmat: struct
%        One capture width matrix struct per computed statistic, keyed by
%        statistic name, e.g. cwmat.mean
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
%     >> [cwmat, maep_matrix] = power_performance_workflow(S, h, P, "mean");
%     maep_matrix =
%         456328960.94  % [W*h]
%
%     >> cwmat.mean
%     ans = 
%       struct with fields:
%            values: [8×14 double]
%              stat: 'mean'
%            x_bins: [0.5000 1.5000 2.5000 3.5000 4.5000 5.5000 … ] (1×14 double)
%            y_bins: [0.2500 0.7500 1.2500 1.7500 2.2500 2.7500 3.2500 3.7500]
%           x_edges: [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]
%           y_edges: [0 0.5000 1 1.5000 2 2.5000 3 3.5000 4]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S
    h
    P
    statistic (1,:) string {mustBeMember(statistic, ["mean", "std", "median", "count", "sum", "min", "max", "probability", "frequency"])}
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

% Bin centers with IEC TS 62600-100 widths of 0.5 m for Hm0 and 1 s for Te,
% so the edges fall on 0, 0.5, 1, ... m and 0, 1, 2, ... s.
Hm0_bins = 0.25:0.5:ceil(max(Hm0) / 0.5) * 0.5;
Te_bins = 0.5:1:ceil(max(Te));

% mean and probability are always needed for MAEP. Any other requested
% statistic is computed once and stored under its own name.
cwmat = struct();
for stat = unique(["mean", "probability", statistic])
    cwmat.(stat) = capture_width_matrix(Hm0, Te, CW, stat, Hm0_bins, Te_bins);
end

jmat = wave_energy_flux_matrix(Hm0, Te, J, "mean", Hm0_bins, Te_bins);

maep_matrix = mean_annual_energy_production_matrix(cwmat.mean, jmat, cwmat.probability);

% Plot each requested statistic
cw_matrix = gobjects(1, numel(statistic));
for i = 1:numel(statistic)
    figure('Name', sprintf('Capture Width Matrix %s', statistic(i)), 'NumberTitle', 'off')
    cw_matrix(i) = plot_matrix(cwmat.(statistic(i)), "Capture Width");
    if strlength(options.savepath) > 0
        name = fullfile(options.savepath, sprintf('Capture Width Matrix %s.png', statistic(i)));
        saveas(cw_matrix(i), name);
    end
end

end
