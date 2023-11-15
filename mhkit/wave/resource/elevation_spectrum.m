function S=elevation_spectrum(ts,sample_rate,nnft,time,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave spectra from wave probe timeseries
%    
% Parameters
% ------------
%     ts: matrix or table
%           Wave probe time-series data, with each column a different time
%           series
%
%     sample_rate: float
%           Data frequency (Hz)
%
%     nnft: integer
%           Number of bins in the Fast Fourier Transform
%
%     time: vector or table
%           time (s)
%
%     window: string (Optional)
%        Signal window type. "hamming" is used by default given the broadband 
%        nature of waves. See scipy.signal.get_window for more options.
%        to call: elevation_spectrum(ts,sample_rate,nnft,time,"window",window)
%
%     detrend: logical (Optional)
%        Specifies if a linear trend is removed from the data before calculating 
%        the wave energy spectrum.  Data is detrended by default.
%        to call: elevation_spectrum(ts,sample_rate,nnft,time,"detrend",detrend)
%
%     noverlap: integer (Optional)
%       Number of points to overlap between segments. If None,
%       noverlap = nperseg / 2.  Defaults to None.
%       to call: elevation_spectrum(ts,sample_rate,nnft,time,"noverlap",noverlap)
%
% Returns
% ---------
%     S: structure
%
%
%         S.spectrum: vector or matrix Spectral Density (m^2/Hz) per probe           
%
%         S.type: 'Spectra from Time Series'
%
%         S.frequency: frequency [Hz] 
%
%         S.sample_rate: sample_rate
%
%         S.nnft: nnft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    ts 
    sample_rate
    nnft
    time
    options.window = "hamming";
    options.detrend = true;
    options.noverlap = py.None;
end

py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');
py.importlib.import_module('pandas');
py.importlib.import_module('mhkit_python_utils');
if (isa(ts,'py.pandas.core.frame.DataFrame')~=1)
    if (isa(ts,'table')==1)
        ts=table2array(ts);
    end
    if (isa(time,'table')==1)
        time=table2array(time);
    end
    x=size(ts);
    li=py.list();
    if x(2)>1 
        for i = 1:x(2)
            app=py.list(ts(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
        end  
        ts=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,time,int32(x(2)));        

    elseif x(2)==1
       ts=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(ts,time,int32(x(2))); 
       
    end
    
end
nnft=int32(nnft);


spectra=py.mhkit.wave.resource.elevation_spectrum(ts,sample_rate,nnft,pyargs('window',options.window,'detrend',options.detrend,'noverlap',options.noverlap));
 vals=double(py.array.array('d',py.numpy.nditer(spectra.values)));
 sha=cell(spectra.values.shape);
 x=int64(sha{1,1});
 y=int64(sha{1,2});
 vals=reshape(vals,[x,y]);


S.spectrum=vals;

S.type='Spectra from Timeseries';

S.frequency=double(py.array.array('d',py.numpy.nditer(spectra.index))).';
S.sample_rate=sample_rate;
S.nnft=nnft;

    
