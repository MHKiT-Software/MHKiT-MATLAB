function out = band_aggregate(spsdl, octave, fmin, fmax, method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reorganize spectral density level frequency tensor into fractional octave bands and applies a function to them
%
% Parameters
% ------------
%   spsdl: struct
%       Mean square sound pressure spectral density level in dB rel 1 uPa^2/Hz
%       spsdl.data : Spectral density level data [dB rel 1 uPa^2/Hz]
%       spsdl.freq : Frequency vector [Hz]
%       spsdl.time : Time vector
%       spsdl.fs : Sampling frequency [Hz]
%   octave: (1,2) vector
%       Octave and octave base to subdivide spectral density level by
%       Set octave base to 2 for true octave band; set to base 10 for decidecade octave band
%       Default = [3, 2] (true third octave)
%   fmin: float
%       Lower frequency band limit (lower limit of the hydrophone) [Hz]
%   fmax: float
%       Upper frequency band limit (Nyquist frequency) [Hz]
%   method: string or cell
%       Method to run on the binned data. Can be a string (e.g., "median") or a cell
%       where the first element is the method and the second is its argument
%       Options: [median, mean, min, max, sum, quantile, std, var, count, map]
%
% Returns
% ---------
%   out: struct
%       Frequency band-averaged sound pressure spectral density level
%       out.data : Band-averaged spectral density level [dB re 1 uPa^2/Hz]
%       out.freq : Center frequencies of bands [Hz]
%       out.time : Time vector
%       out.name : Data name
%       out.units : Data units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% Validate spsdl structure
validate_spsdl_struct(spsdl, 'band_aggregate', ...
    'required_fields', {{'data', 'freq', 'time', 'name', 'units', 'fs'}});

if ~isnumeric(octave) || numel(octave) ~= 2
    error('MHKiT:acoustics:band_aggregate:InvalidOctave', 'octave must be a vector of two positive integers');
end
if any(octave <= 0)
    error('MHKiT:acoustics:band_aggregate:InvalidOctave', 'octave must contain positive integers');
end
if ~isnumeric(fmin) || fmin <= 0
    error('MHKiT:acoustics:band_aggregate:InvalidFrequency', 'fmin must be a positive number');
end
if ~isnumeric(fmax) || fmax <= fmin
    error('MHKiT:acoustics:band_aggregate:InvalidFrequency', 'fmax must be greater than fmin');
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
