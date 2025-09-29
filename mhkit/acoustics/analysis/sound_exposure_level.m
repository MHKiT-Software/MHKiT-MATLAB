function sel = sound_exposure_level(spsd, group, fmin, fmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate sound exposure level (SEL) across a specified frequency band
% from the sound pressure spectral density (SPSD). If a marine mammal group is
% provided, the resulting SEL is weighted according to the U.S. National Marine
% Fisheries Service (NMFS) guidelines.
%
% Parameters
% ------------
%   spsd: struct
%       spsd.data : Sound pressure spectral density [Pa^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector [s]
%       spsd.fs : Sampling frequency [Hz]
%       spsd.nfft : Number of FFT points
%       spsd.nbin : Bin length for exposure computation
%   group: char or string
%       Marine mammal group for auditory weighting function. Options: 
%         'LF' (low frequency cetaceans),
%         'HF' (high frequency cetaceans), 
%         'VHF' (very high frequency cetaceans),
%         'PW' (phocid pinnepeds), 
%         'OW' (otariid pinnepeds). 
%         Default: []
%   fmin: double
%       Lower frequency band limit [Hz]. Default: 10
%   fmax: double
%       Upper frequency band limit [Hz]. Default: 100000
%
% Returns
% ---------
%   sel: struct
%       sel.data : Sound exposure level [dB re 1 uPa^2 s]
%       sel.time : Time vector [s]
%       sel.units : Units string
%       sel.name : Descriptive name
%       sel.group : Marine mammal group used
%       sel.integration_time : Integration time [s]
%       sel.freq_band_min : Lower frequency limit [Hz]
%       sel.freq_band_max : Upper frequency limit [Hz]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    group = []
    fmin double {mustBeFinite, mustBePositive} = 10
    fmax double {mustBeFinite, mustBePositive} = 100000
end

arguments (Output)
    sel struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'sound_exposure_level', ...
    'required_fields', {{'freq', 'time', 'fs', 'nfft', 'nbin'}});

% Validate frequency bounds
if fmax <= fmin
    error('MHKiT:acoustics:sound_exposure_level:InvalidInput', 'fmax must be greater than fmin');
end

% Validate group parameter if provided
if ~isempty(group)
    valid_groups = {'LF', 'HF', 'VHF', 'PW', 'OW'};
    if ~ismember(group, valid_groups)
        error('MHKiT:acoustics:sound_exposure_level:InvalidInput', ...
            'group must be one of: %s', strjoin(valid_groups, ', '));
    end
end
fmax = fmax_warning(spsd.fs/2, fmax);

% get weighting factor
if ~isempty(group)
    [w, e] = nmfs_auditory_weighting(spsd.freq, group);
    w = 10 .^ (w/10);
    name = "Weighted Sound Exposure Level";
else
    w = ones(length(spsd.freq), 1);
    name = "Sound Exposure Level";
end

% Reference value of sound pressure
reference = 1e-12 * 1;  % Pa^2 s, = 1 uPa^2 s

% Select frequency band
idx = spsd.freq >= fmin & spsd.freq <= fmax;
band = spsd.data(idx, :);
freqs = spsd.freq(idx);
w = w(idx);

exposure = trapz(freqs, band.*w);

% Sound exposure level (L_{E,p}) = (L_{p,rms} + 10log10(t))
s = 10 * log10(exposure / reference) + 10 * log10(spsd.nfft / spsd.fs);  % n_points / (n_points/s)

sel.data = s;
sel.time = spsd.time;
sel.units = "dB re 1 uPa^2 s";
sel.name = name;
sel.group = group;
sel.integration_time = spsd.nbin;
sel.freq_band_min = fmin;
sel.freq_band_max = fmax;


end
