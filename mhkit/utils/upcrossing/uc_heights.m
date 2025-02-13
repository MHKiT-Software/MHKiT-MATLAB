function heights = uc_heights(t, data, inds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the height between zero crossings.
%
% The height is defined as the max value - min value
% between the zero crossing points.
%
% Parameters
% ------------
%     t: array
%         Array of time values.
%     data: array
%         Array of data values.
%     inds: array, optional
%         Indices for the upcrossing.
%
% Returns:
% ------------
% heights: Height values of the time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    inds = upcrossing(t, data);
end

heights = uc_apply_(t, data, @(ind1, ind2) max(data(ind1:ind2)) - min(data(ind1:ind2)), inds);
end

