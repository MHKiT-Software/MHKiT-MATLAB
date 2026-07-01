function [weighting_func,exposure_func] = nmfs_auditory_weighting(frequency, group)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the auditory weighting and exposure functions for marine mammals
% based on the National Marine Fisheries Service (NMFS) guidelines.
% 
% The weighting function is applied to sound exposure level to determine the
% auditory impact on marine mammals. The exposure function is the inverse of the
% weighting function and illustrates how the weighting function relates to marine
% mammal hearing thresholds.
% Both function are returned in their log10-transform, in units of dB. To transform
% back to linear units, use 10**(weighting_func/10).
% 
% https://www.fisheries.noaa.gov/national/marine-mammal-protection/marine-mammal-acoustic-technical-guidance-other-acoustic-tools
% 
% Parameters
% ----------
% frequency: array
%     Frequency vector in [Hz].
% group: str
%     Marine mammal group for which the auditory weighting function is applied.
%     Options: 'LF' (low frequency cetaceans), 'HF' (high frequency cetaceans),
%     'VHF' (very high frequency cetaceans), 'PW' (phocid pinnepeds),
%     'OW' (otariid pinnepeds)
% 
% Returns
% -------
% weighting_func: struct
%     Auditory weighting function [unitless] indexed by frequency
% exposure_func: struct
%     Log-transformed auditory exposure function [dB] indexed by frequency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    frequency {mustBeVector}
    group {mustBeText}
end

arguments (Output)
    weighting_func {mustBeVector}
    exposure_func {mustBeVector}
end

group = lower(group);
valid_groups = {'lf', 'hf', 'vhf', 'pw', 'ow'};
if ~ismember(group, valid_groups)
    error('Group must be one of: LF, HF, VHF, PW, OW');
end

params = struct( ...
    'lf', struct('a', 0.99, 'b', 5, 'f1', 0.168, 'f2', 26.6, 'c', 0.12, 'k', 177), ...
    'hf', struct('a', 1.55, 'b', 5, 'f1', 1.73, 'f2', 129, 'c', 0.32, 'k', 181), ...
    'vhf', struct('a', 2.23, 'b', 5, 'f1', 5.93, 'f2', 186, 'c', 0.91, 'k', 160), ...
    'pw', struct('a', 1.63, 'b', 5, 'f1', 0.81, 'f2', 68.3, 'c', 0.29, 'k', 175), ...
    'ow', struct('a', 1.58, 'b', 5, 'f1', 2.53, 'f2', 43.8, 'c', 1.37, 'k', 178) ...
);

p = params.(group);
frequency_kHz = frequency / 1000; % Convert to kHz
ratio_a = frequency_kHz / p.f1;
ratio_b = frequency_kHz / p.f2;
band_filter = ratio_a .^ (2 * p.a) ./ ( ((1 + ratio_a.^2) .^ p.a) .* ((1 + ratio_b.^2) .^ p.b) );

weighting_func = p.c + 10 * log10(band_filter); % dB
exposure_func = p.k - 10 * log10(band_filter); % dB


end