function out = time_aggregate(spsdl, window, method)
% 
% Reorganizes spectral density level frequency tensor into
% time windows and applies a function to them.
% 
% Parameters
% ----------
% spsdl: struct
%     Mean square sound pressure spectral density level in dB rel 1 uPa^2/Hz
% window: int
%     Time in seconds to subdivide spectral density level into. Default: 60 s.
% method: 
%     Method to run on the binned data. Can be a string (e.g., "median") or a dict
%     where the key is the method and the value is its argument (e.g., {"quantile": 0.25}).
%     Options: [median, mean, min, max, sum, quantile, std, var, count]
% 
% Returns
% -------
% out: xarray.DataArray (time_bins, freq)
%     Time-averaged sound pressure spectral density level [dB re 1 uPa^2/Hz]
%     indexed by time and frequency

arguments (Input)
    spsdl struct
    window {mustBeInteger} = 60 %s
    method = "median"
end

arguments (Output)
    out struct
end

if window <= 0
    error("'window' must be a positive integer.")
end

% validate method
[method_name, method_arg] = validate_method(method);
if isempty(method_arg)
    mfunc = str2func(method_name);
elseif method_arg == "quantile"
    mfunc = @(x)quantile(x,method_arg);
else
    mfunc = method_arg;
end

% groupby and apply method
tbins = spsdl.time(1):seconds(window):spsdl.time(end);
idx_bin = discretize(spsdl.time, tbins);
temp = zeros(length(spsdl.freq),length(tbins)-1);
for x=1:length(spsdl.freq)
    temp(x,:) = splitapply(mfunc,spsdl.data(x,:),idx_bin);
end
out.data = temp;
out.time = tbins(1:end-1) + seconds(window)/2;
out.freq = spsdl.freq;
out.units = spsdl.units;
out.name = spsdl.name;

end