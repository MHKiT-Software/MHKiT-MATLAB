function [t_peaks,peaks] = global_peaks(t,data)

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

if ~isa(t,'numeric')
    error('ERROR: t must be a double array')
end
if ~isa(data,'numeric')
    error('ERROR: data must be a double array')
end

py.importlib.import_module('mhkit');

t = py.numpy.array(t);
data = py.numpy.array(data);

result = py.mhkit.loads.extreme.global_peaks(t, data);

t_peaks = double(result{1});
peaks = double(result{2});

end



