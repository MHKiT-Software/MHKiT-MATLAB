function troughs = uc_troughs(t, data, inds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the troughs between zero crossings.
%
% Parameters
% ----------
% t : array
%     Time array.
% data : array
%     Signal time-series.
% inds : array, optional
%     Indices for the upcrossing.
%
% Returns
% -------
% troughs : array
%     Trough values of the time-series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    inds = upcrossing(t, data);
end

troughs = uc_apply_(t, data, @(ind1, ind2) min(data(ind1:ind2)), inds);
end


