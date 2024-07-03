function block_maxima = block_maxima(t, x, t_st)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Find the block maxima of a time-series.
%
%     The timeseries (t,x) is divided into blocks of length t_st, and the 
%     maxima of each block is returned.
%
%     Parameters
%     ----------
%         t : array
%             Time array
%         x : array
%             Global peaks time-series
%         t_st : double
%             Short-term period
%
%      Returns
%      -------
%         block_maxima : array
%             Block maxima (i.e. largest peak in each block)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(t,'numeric')
    error('ERROR: t must be a double array')
end
if ~isa(x,'numeric')
    error('ERROR: x must be a double array')
end
if ~isa(t_st,'double')
    error('ERROR: t_st must be a double')
end

py.importlib.import_module('mhkit');

t = py.numpy.array(t);
x = py.numpy.array(x);

result = py.mhkit.loads.extreme.block_maxima(t, x, t_st);

block_maxima = double(result);

end