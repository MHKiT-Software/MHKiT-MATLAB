function out = time_aggregate(spsdl, window, method)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reorganize spectral density level frequency tensor into time windows and apply a function to them.
%
% Parameters
% ------------
%   spsdl: structure
%       spsdl.data : Sound pressure spectral density level data [dB re 1 uPa^2/Hz]
%       spsdl.time : Time vector
%       spsdl.freq : Frequency vector [Hz]
%       spsdl.units : Data units string
%       spsdl.name : Data name string
%   window: integer
%       Time in seconds to subdivide spectral density level into. Default: 60 s.
%   method: string or structure
%       Method to run on the binned data. Can be a string (e.g., "median") or a structure
%       where the field is the method and the value is its argument (e.g., struct('quantile', 0.25)).
%       Options: [median, mean, min, max, sum, quantile, std, var, count]
%
% Returns
% ---------
%   out: structure
%       out.data : Time-averaged sound pressure spectral density level [dB re 1 uPa^2/Hz]
%       out.time : Time vector for binned data
%       out.freq : Frequency vector [Hz]
%       out.units : Data units string
%       out.name : Data name string
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsdl struct
    window {mustBeInteger} = 60 % seconds
    method = "median"
end

arguments (Output)
    out struct
end

% Validate spsdl structure
validate_spsdl_struct(spsdl, 'time_aggregate');

if window <= 0
    error('MHKiT:acoustics:time_aggregate:InvalidInput', ...
        'window must be a positive integer');
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
