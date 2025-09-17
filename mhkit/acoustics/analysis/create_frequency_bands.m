function [octave_bins, band] = create_frequency_bands(octave, base, fmin, fmax)

% Calculates frequency bands based on the specified octave, minimum and
% maximum frequency limits.
%
% Parameters
% ----------
% octave: int
%     Octave to subdivide spectral density level by.
% base : int, optional
%     Octave base. Set to 2 for the true octave band; set to base 10 for
%     the decidecade octave band. Default: 2
% fmin : int, optional
%     Lower frequency band limit (lower limit of the hydrophone). Default is 10 Hz.
% fmax : int, optional
%     Upper frequency band limit (Nyquist frequency). Default is 100,000 Hz.
%
% Returns
% -------
% octave_bins: array
%     Array of octave bin edges
% band: struct
%     Struct containing the frequency band edges and center frequency

bandwidth = base^(1 / octave);
half_bandwidth = base^(1 / (octave * 2));

% Calculate center frequencies
log_fmin = log10(fmin);
log_fmax = log10(fmax * bandwidth);
step = log10(bandwidth);
center_freq = 10 .^ (log_fmin : step : log_fmax);

lower_limit = center_freq / half_bandwidth;
upper_limit = center_freq * half_bandwidth;
octave_bins = [lower_limit, upper_limit(end)];

band = struct();
band.center_freq = center_freq;
band.lower_limit = lower_limit;
band.upper_limit = upper_limit;

end