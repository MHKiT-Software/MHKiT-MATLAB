function sel = sound_exposure_level(spsd, group, fmin, fmax)
% 
% Calculates the sound exposure level (SEL) across a specified frequency band
% from the sound pressure spectral density (SPSD). If a marine mammal group is
% provided, the resulting SEL is weighted according to the U.S. National Marine
% Fisheries Service (NMFS) guidelines.
% 
% Parameters
% ----------
% spsd: struct
%     Sound pressure spectral density in [Pa^2/Hz] with a bin length
%     equal to the time over which sound exposure should be computed.
% group: str
%     Marine mammal group for which the auditory weighting function is applied.
%     Options: 'LF' (low frequency cetaceans), 'HF' (high frequency cetaceans),
%     'VHF' (very high frequency cetaceans), 'PW' (phocid pinnepeds),
%     'OW' (otariid pinnepeds). Default: None
% fmin: int
%     Lower frequency band limit (lower limit of the hydrophone).
%     Default: 10 Hz
% fmax: int
%     Upper frequency band limit (Nyquist frequency). Default:
%     100000 Hz
% 
% Returns
% -------
% sel: struct
%     Sound exposure level [dB re 1 uPa^2 s] indexed by time

arguments (Input)
    spsd struct
    group = []
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
end

arguments (Output)
    sel struct
end

% type checks
if fmin <= 0
    error('fmin must be positive');
end
if fmax <= fmin
    error('fmax must be greater than fmin');
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