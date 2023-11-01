function out = block_maxima(t, x, t_st)
%     Find the block maxima of a time-series.
% 
%     The timeseries (t,x) is divided into blocks of length t_st, and the
%     maxima of each bloock is returned.
% 
%     Parameters
%     ----------
%     t : array
%         Time array.
%     x : array
%         global peaks timeseries.
%     t_st : float
%         Short-term period.
% 
%     Returns
%     -------
%     block_maxima: array
%         Block maxima (i.e. largest peak in each block).

assert(isvector(t),'t must be array')
assert(isvector(x),'data must be array')
assert(isfloat(t_st),'t_st must be of type float')

nblock = int32(t(end)/t_st);
out = zeros([1,nblock]);
for i = 1:nblock
    ix = x((t >= (i-1)*t_st) & (t < i*t_st));
    out(i) = max(ix);
end

end

