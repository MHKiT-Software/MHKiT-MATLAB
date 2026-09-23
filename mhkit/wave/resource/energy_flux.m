function J = energy_flux(S, h, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the omnidirectional wave energy flux of the spectra
%
% For deep water (deep=true), uses simplified formula (Eq 8 in
% IEC 62600-100 Ed. 2.0 en 2024):
%   J = (rho * g^2 / (64*pi)) * Hm0^2 * Te
%
% For general depth (deep=false), uses Eq 9 in IEC 62600-101 Ed. 2.0 en 2024:
%   J = rho * g * sum(Cg * S * df)
%
% Parameters
% ------------
% S : struct
%   Wave spectrum structure:
%     S.spectrum : vector or matrix [m^2/Hz]
%       Spectral density
%     S.frequency : vector [Hz]
%       Frequency
%     S.time : datetime or duration (optional)
%       One timestamp per spectrum/column
%   Numeric, table, and timetable input are not supported here: h (water
%   depth) occupies the positional slot the multi-container convention
%   uses for an explicit frequency argument, so only struct input (which
%   carries frequency inline) is unambiguous.
% h : double [m]
%   Water depth
% deep : logical (optional)
%   If true, use the deep water approximation. Default = false.
%   When false, a depth check is run to check for shallow water.
% rho : double [kg/m^3] (optional)
%   Water density. Default = 1025 kg/m^3
% g : double [m/s^2] (optional)
%   Gravitational acceleration. Default = 9.80665 m/s^2
% ratio : double (optional)
%   Only applied if deep=false. If h/l > ratio,
%   water depth is set to deep. Default = 2
%
% Returns
% ---------
% J : double or column vector [W/m]
%   Omni-directional wave energy flux, one value per spectrum. A table
%   or timetable if S.time is present, matching mhkit_restore_spectrum_output.
%
% Examples
% --------
%     CDIP Example: Station Number 225: Kaneohe Bay, WETS, Oahu, HI
%     https://cdip.ucsd.edu/themes/cdip?pb=1&u2=s:225:st:1&d2=p70
%     >> station_number = '225';
%     >> data_type = 'realtime';
%     >> years = 2025;
%     >> parameters = {'waveEnergyDensity', 'metaWaterDepth'};
%     >> data = cdip_request_parse_workflow('station_number', station_number, 'data_type', data_type, 'years', years, 'parameters', parameters);
%     >> % cast to double: fzero (inside wave_number, called by energy_flux
%     >> % below) requires double, but CDIP returns these as single.
%     >> frequency = double(data.metadata.wave.waveFrequency);  % [Hz], 64x1 double
%     >> % waveEnergyDensity is stored [time x frequency]; transpose to
%     >> % match MHKiT's frequency-as-rows, spectra-as-columns convention.
%     >> spectrum = data.data.wave2D.waveEnergyDensity';  % [m^2/Hz], 64x17520 double
%     >> % CDIP's own buoy depth reading, not a literature value.
%     >> h = double(data.metadata.meta.metaWaterDepth);  % [m], 84
%
%     Struct (CDIP real-world data): single (most recent) spectrum
%     >> S.frequency = frequency;
%     >> S.spectrum = spectrum(:,end);
%     >> J = energy_flux(S, h);
%     J =
%         14495.52  % [W/m]
%
%     Struct (CDIP real-world data): full dataset, one value per spectrum
%     >> S.spectrum = spectrum;
%     >> J = energy_flux(S, h);
%     J =  % [W/m]
%       10152.76
%        9548.60
%        9148.13
%          :
%       16898.06
%       20193.00
%       14495.52
%
%     WEC-Sim Output Example
%     >> S = load('examples/data/RM3MooringMatrix_matlabWorkspace.mat', 'output');
%     >> elevation = S.output.wave.elevation;  % [m], RM3 float, 40001x1 double
%     >> raw_time = S.output.wave.time;  % [s], 40001x1 double
%     >> sample_rate = 1 / (raw_time(2) - raw_time(1));  % [Hz], 100
%     >> Sxx = elevation_spectrum(elevation, sample_rate, 1000, raw_time);
%     >> frequency = Sxx.frequency;  % [Hz], 501x1 double
%     frequency =
%         0.0000
%         0.1000
%         0.2000
%          :
%        49.8000
%        49.9000
%        50.0000
%
%     >> spectrum = Sxx.spectrum;  % [m^2/Hz], 501x1 double
%     spectrum =
%         0.1439
%         0.8485
%         0.8739
%          :
%         0.0000
%         0.0000
%         0.0000
%
%     >> % The official WEC-Sim RM3 tutorial's hydrodynamic data was
%     >> % computed by WAMIT at infinite depth (rm3.out: "Water depth:
%     >> % infinite"), so this example uses the deep-water approximation.
%     >> % h is unused when deep=true; Inf documents that assumption
%     >> % rather than standing in for a real finite depth.
%     Struct (WEC-Sim output): the single spectrum, deep-water approximation
%     >> S.frequency = frequency;
%     >> S.spectrum = spectrum;
%     >> J = energy_flux(S, Inf, 'deep', true);
%     J =
%         10741.91  % [W/m]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    S struct
    h (1,1) {mustBeNumeric}
    options.deep (1,1) logical = false
    options.rho (1,1) {mustBeNumeric, mustBePositive} = 1025
    options.g (1,1) {mustBeNumeric, mustBePositive} = 9.80665
    options.ratio (1,1) {mustBeNumeric, mustBePositive} = 2
end

arguments (Output)
    J
end

[spectrum, frequency, time, input_style] = mhkit_standardize_spectrum_input(S, 'energy_flux');

rho = options.rho;
g = options.g;

% Validate depth is positive when using general depth calculation
if ~options.deep && h <= 0
    error('MHKiT:energy_flux:InvalidInput', ...
        'h must be positive when deep=false');
end

if options.deep
    % Eq 8 in IEC 62600-100 Ed. 2.0 en 2024 (deep water simplification)
    Te = energy_period(spectrum, frequency);
    Hm0 = significant_wave_height(spectrum, frequency);

    coeff = rho * (g^2) / (64 * pi);
    J = coeff * (Hm0.^2) .* Te;
else
    % Calculate wave number
    k = wave_number(frequency, h, 'rho', rho, 'g', g);

    % Calculate wave celerity (group velocity)
    Cg = wave_celerity(k, h, 'g', g, 'depth_check', true, 'ratio', options.ratio);

    % Calculate frequency bin widths
    delta_f = diff(frequency);
    delta_f = [frequency(2) - frequency(1); delta_f];  % Prepend first bin width

    % Eq 9 in IEC 62600-101 Ed. 2.0 en 2024
    Cg_values = Cg.values(:);

    if isvector(spectrum)
        J = rho * g * sum(Cg_values .* spectrum .* delta_f);
    else
        % Matrix case: sum along frequency dimension (first dimension)
        J = rho * g * sum(Cg_values .* spectrum .* delta_f, 1);
    end
end

J = J(:);
mhkit_verify_is_column_vector(J, 'energy_flux');
J = mhkit_restore_spectrum_output(J, input_style, 'energy_flux', time);

end
