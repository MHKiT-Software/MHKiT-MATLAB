function wave_elevation=surface_elevation(S,time_index,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave elevation time series from spectrum using a random phase
%    
% Parameters
% ------------
%    S: Spectral Density (m^2/Hz)
%       structure of form:
%           S.spectrum: Spectral Density (m^2/Hz)
%
%           S.type: String of the spectra type, i.e. Bretschneider, 
%           time series, date stamp etc.
%
%           S.frequency: frequency (Hz)
%
%    time_index: array
%        Time used to create the wave elevation time series [s]
%        
%    seed: Int (optional)
%        random seed
%        to call: wave_elevation(S,time_index,"seed",seed)
%
%    frequency_bins: vector (optional) 
%       Bin widths for frequency of S. Required for unevenly sized bins
%       to call: wave_elevation(S,time_index,"frequency_bins",frequency_bins)
%
%    phases: vector or matrix (optional)
%       Explicit phases for frequency components (overrides seed)
%       to call: wave_elevation(S,time_index,"phases",phases)
%     
% Returns
% ---------
%    wave_elevation: structure
%
%
%         wave_elevation.elevation: Wave surface elevation (m)
%
%         wave_elevation.type: 'Time Series from Spectra'
%
%         wave_elevation.time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    S 
    time_index
    options.seed {mustBeNumeric} = 123;
    options.frequency_bins = nan;
    options.phases = nan;
end

if ~isnan(options.frequency_bins)
    assert(all(size(options.frequency_bins) == size(S.frequency)), ...
        "shape of frequency_bins must match shape of S")
end

if ~isnan(options.phases)
    assert(all(size(options.phases) == size(S.frequency)), "shape " + ...
        "of phases must match shape of S")
end


f = S.frequency;

if isnan(options.frequency_bins)
    delta_f = f(2)-f(1);
else
    delta_f = options.frequency_bins;
end

if isnan(options.phases)
    rng(options.seed);
    phase = 2*pi*rand(size(S.spectrum));
else
    phase = options.phases;
end

omega = 2*pi*f;

% Wave amplitude times delta f
A = 2*S.spectrum;
A = A.*delta_f;
A = sqrt(A);

% Product of omega and time
B = time_index.*omega;

% wave elevation
C = cos(B + phase);
elevation = sum((C.*A));

wave_elevation.type = 'Time Series from Spectra';
wave_elevation.elevation = elevation;
wave_elevation.time = time_index;