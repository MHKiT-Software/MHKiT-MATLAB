function peaks = uc_peaks(t, data, inds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds the peaks between zero crossings.
%
% Parameters:
% ------------
%   t: array
%       Time array.
%   data: array
%       Signal time-series.
%   inds: Optional, array
%       indices for the upcrossing.
%
% Returns:
% ------------
%   peaks: array
%       Peak values of the time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    inds = upcrossing(t, data);
end

peaks = uc_apply_(t, data, @(ind1, ind2) max(data(ind1:ind2)), inds);
end

