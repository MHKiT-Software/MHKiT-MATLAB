function frequency=instantaneous_frequency(voltage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calcultes the instantaneous frequency of a measured voltage

% Parameters
% -----------
%     voltage: structure with handles- voltage.voltage and voltage.time  
%         measured voltage source (V) index by time 
% 
% Returns
% -------
%     THCD: double
%         the total harmonic current distortion
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

time= voltage.time ;
data = voltage.voltage;
dname = 'voltage';

dsize=size(data);

li=py.list();
if dsize(2)>1 
   for i = 1:dsize(2)
      app=py.list(data(:,i));
      li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
   end
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time(:,1),li,int32(dsize(2)));
elseif dsize(2)==1
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,py.numpy.array(data),dsize(x(2)));
end

harmonics_pd = py.mhkit.power.characteristics.instantaneous_frequency(data_pd);

vals=double(py.array.array('d',py.numpy.nditer(harmonics_pd.values)));
sha=cell(harmonics_pd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[x,y]);

frequency.frequency = vals;
frequency.time = time; 

