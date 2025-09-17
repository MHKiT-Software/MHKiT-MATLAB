function out = band_aggregate(spsdl, octave, fmin, fmax, method)

% Reorganizes spectral density level frequency tensor into
% fractional octave bands and applies a function to them.
% 
% Parameters
% ----------
% spsdl: struct
%     Mean square sound pressure spectral density level in dB rel 1 uPa^2/Hz
% octave: [int, int]
%     Octave and octave base to subdivide spectral density level by. Set to
%     octave base to 2 for the true octave band; set to base 10 for
%     the decidecade octave band.
%     Default = [3, 2] (true third octave)
% fmin: int
%     Lower frequency band limit (lower limit of the hydrophone). Default: 10 Hz
% fmax: int
%     Upper frequency band limit (Nyquist frequency). Default: 100000 Hz
% method: 
%     Method to run on the binned data. Can be a string (e.g., "median") or a dict
%     where the key is the method and the value is its argument (e.g., {"quantile": 0.25}).
%     Options: [median, mean, min, max, sum, quantile, std, var, count, map]
% 
% Returns
% -------
% out: struct
%     Frequency band-averaged sound pressure spectral density level [dB re 1 uPa^2/Hz]
%     indexed by time and frequency


arguments (Input)
    spsdl struct
    octave {mustBeVector} = [3,2]
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
    method = "median"
end

arguments (Output)
    out struct
end

if ~isnumeric(octave) || numel(octave) ~= 2
    error("'octave' must be a vector of two positive integers.");
end
if any(octave <= 0)
    error("'octave' must contain positive integers.");
end
if ~isnumeric(fmin) || fmin <= 0
    error("'fmin' must be a positive integer.");
end
if ~isnumeric(fmax) || fmax <= fmin
    error("'fmax' must be greater than 'fmin'.");
end

fmax = fmax_warning(spsdl.fs/2, fmax);

% validate method
[method_name, method_arg] = validate_method(method);
if isempty(method_arg)
    mfunc = str2func(method_name);
elseif method_arg == "quantile"
    mfunc = @(x)quantile(x,method_arg);
else
    mfunc = method_arg;
end

% derive octave bins
[octave_bins, band] = create_frequency_bands(octave(1),octave(2),fmin,fmax);

% groupby and apply method
idx_bin = discretize(spsdl.freq, octave_bins);
temp = zeros(length(band.center_freq), length(spsdl.time));
for x=1:length(spsdl.time)
    temp(:,x) = splitapply(mfunc,spsdl.data(:,x),idx_bin);
end

out.data = temp;
out.freq = band.center_freq(min(idx_bin):max(idx_bin));
out.time = spsdl.time;
out.name = spsdl.name;
out.units = spsdl.units;

end