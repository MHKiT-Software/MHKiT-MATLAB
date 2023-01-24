function out = fill_time_gaps(epoch, sample_rate_hz)
% Fill gaps (NaN values) in the timeseries by simple linear
% interpolation.  The ends are extrapolated by stepping
% forward/backward by 1/sample_rate_hz.

dt = 1 / sample_rate_hz;
% using fillgaps to interplolate values which uses nan so we convert 
% negative values to nan
missing_data = epoch < 0;
epoch(missing_data) = nan;
out = fillgaps(epoch, false);
if isnan(epoch(1))
    i0 = find(~isnan(epoch),1) - 1;
    delta = (-i0:-1)*dt;
    epoch(1:i0) = epoch(i0) + delta;
end
if isnan(epoch(end))
    ie = find(~isnan(epoch),1,'last');
    delta = (1:numel(epoch)-ie) * dt;
    epoch(ie+1:end) = epoch(ie) + delta;
end

end

