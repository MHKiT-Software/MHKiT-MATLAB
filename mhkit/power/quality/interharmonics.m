function interharmonics=interharmonics(harmonics,grid_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the interharmonics from the harmonics of current based on IEC 61000-4-7.
%
% Parameters
% -----------
%     harmonics: structure with handles- harmonics.amplitude and harmonics.harmonic
%         Harmonic amplitude with each timeseries in its own column
%
%     grid_freq: int
%         Value indicating if the power supply is 50 or 60 Hz. Options = 50 or 60
%
% Returns
% -------
%     interharmonics: structure with handles interharmonics.amplitude and
%           interharmonics.harmonic
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

data = harmonics.amplitude;

dsize=size(data);
li=py.list();
if dsize(2)>1
   for i = 1:dsize(2)
      app=py.list(data(:,i));
      li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

   end
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(harmonics.harmonic,li,int32(dsize(2)));
elseif dsize(2)==1
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(harmonics.harmonic,py.numpy.array(data),int32(dsize(2)));
end

interharmonics_pd = py.mhkit.power.quality.interharmonics(data_pd,grid_freq);
vals=double(py.array.array('d',py.numpy.nditer(interharmonics_pd.values)));
sha=cell(interharmonics_pd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[x,y]);


interharmonics.amplitude=vals;
interharmonics.harmonic = double(py.array.array('d',py.numpy.nditer(interharmonics_pd.index)));

