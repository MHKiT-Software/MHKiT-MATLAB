function [octave_bins, band] = create_frequency_bands(octave, base, fmin, fmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate frequency bands based on the specified octave, minimum and maximum frequency limits
%
% Parameters
% ------------
%   octave: double
%       Octave to subdivide spectral density level by
%   base: double, optional
%       Octave base. Set to 2 for the true octave band; set to base 10 for
%       the decidecade octave band. Default: 2
%   fmin: double, optional
%       Lower frequency band limit (lower limit of the hydrophone). Default is 10 Hz
%   fmax: double, optional
%       Upper frequency band limit (Nyquist frequency). Default is 100,000 Hz
%
% Returns
% ---------
%   octave_bins: array
%       Array of octave bin edges
%   band: structure
%       band.center_freq : Center frequencies [Hz]
%       band.lower_limit : Lower frequency limits [Hz]
%       band.upper_limit : Upper frequency limits [Hz]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    octave double {mustBeFinite, mustBePositive}
    base double {mustBeFinite, mustBePositive} = 2
    fmin double {mustBeFinite, mustBePositive} = 10
    fmax double {mustBeFinite, mustBePositive} = 100000
end

arguments (Output)
    octave_bins double
    band struct
end

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
