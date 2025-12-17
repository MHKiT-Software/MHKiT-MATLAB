function [t_peaks,peaks] = global_peaks(time,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Find the global peaks of a zero-cenered response time-series.
%
%     The global peaks are the maxima between consecutive zero
%     up-crossings.
%
%     Parameters
%     ----------
%         t : array
%             Time array
%         data: array
%             Response time-series
%
%      Returns
%      -------
%         t_peaks : array
%             Time array for peaks
%         peaks : array
%             Peak values of the response time-series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input validation
if ~isa(time,'numeric')
    error('ERROR: t must be a double array')
end
if ~isa(data,'numeric')
    error('ERROR: data must be a double array')
end
if length(time) ~= length(data)
    error('time and data must have the same length');
end

% Find zero up-crossings
inds = upcrossing(time, data);

% Include the final point in the dataset
inds = [inds; length(data)];

% Find peak indices between consecutive zero up-crossings
peak_inds = zeros(length(inds)-1, 1);

for i = 1:length(inds)-1
    ind1 = inds(i);
    ind2 = inds(i+1);

    % Find the index of maximum value in this segment
    [~, local_max_idx] = max(data(ind1:ind2-1));
    peak_inds(i) = ind1 + local_max_idx - 1;
end

% Return peak times and values
t_peaks = time(peak_inds);
peaks = data(peak_inds);

end


