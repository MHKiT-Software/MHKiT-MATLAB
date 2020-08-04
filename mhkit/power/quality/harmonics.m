function harmonics=harmonics(x,freq,grid_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the harmonics from time series of voltage or current based on IEC 61000-4-7.
%    
% Parameters
% -----------
%     x: structure with x.time and x.current or x.voltage as values 
%         Time-series of voltage [V] or current [A] 
% 
%     freq: double
%         Frequency of the time-series data [Hz]
%
%     grid_freq: int
%         Value indicating if the power supply is 50 or 60 Hz. Options = 50 or 60
% 
% Returns
% -------
%     harmonics: structure
%         harmonic amplitude and frequency of the time-series data
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

time= x.time ;

if isfield(x, 'current')
    data = x.current;
    dname = 'current';
elseif isfield(x,'voltage')
    data = x.voltage;
    dname = 'voltage';
else 
    ME = MException('MATLAB:harmonics','invalid handles in structure, must contain x.current or x.voltage');
        throw(ME);
end

dsize=size(data);

li=py.list();
if dsize(2)>1 
   for i = 1:dsize(2)
      app=py.list(data(:,i));
      li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
   end
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time(:,1),li,int32(dsize(2)));
elseif dsize(2)==1
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,py.numpy.array(data),dsize(2));
end

harmonics_pd = py.mhkit.power.quality.harmonics(data_pd,freq,grid_freq);

vals=double(py.array.array('d',py.numpy.nditer(harmonics_pd.values)));
sha=cell(harmonics_pd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[x,y]);



harmonics.amplitude=vals;
harmonics.harmonic = double(py.array.array('d',py.numpy.nditer(harmonics_pd.index)));
harmonics.type = dname;


