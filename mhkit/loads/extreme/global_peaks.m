function [t_peaks,peaks] = global_peaks(t,data)
% Find the global peaks of a zero-centered response time-series.
% 
% The global peaks are the maxima between consecutive zero
% up-crossings.
% 
% Parameters
% ----------
% t: array
%     Time array.
% data: array
%     Response time-series.
% 
% Returns
% -------
% t_peaks: array
%     Time array for peaks
% peaks: array
%     Peak values of the response time-series

if ~isvector(t)
    ME = MException('MATLAB:extreme:global_peaks',...
        'Time must be a vector');
    throw(ME);
end
if ~isvector(data)
    ME = MException('MATLAB:extreme:global_peaks',...
        'Data must be a vector');
    throw(ME);
end

% Eliminate Zeros
data(data == 0) = 0.5 * min(abs(data));

% Zero up-crossing
dif = diff(sign(data));
zeroUpCrossings_mask = dif == 2 | dif == 1;
zeroUpCrossings_index = find(zeroUpCrossings_mask);
zeroUpCrossings_index(end+1) = length(data);

% global peaks
npeaks = length(zeroUpCrossings_index);
peaks = zeros([1,npeaks-1]);
t_peaks = zeros([1,npeaks-1]);
for i = 1:npeaks-1
    [peak, peak_index] = max(...
        data(zeroUpCrossings_index(i):zeroUpCrossings_index(i+1)));
    t_peaks(i) = t(zeroUpCrossings_index(i) + peak_index - 1);
    peaks(i) = peak;
end

end

