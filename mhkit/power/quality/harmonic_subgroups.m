function harmonic_subgroups=harmonic_subgroups(harmonics,grid_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the harmonic subgroups from the harmonics of current
%
% Parameters
% -----------
%     harmonics: structure with handles- harmonics.amplitude and harmonics.harmonic
%         harmonic amplitude with each timeseries in its own column
%
%     grid_freq: int
%         value indicating if the power supply is 50 or 60 Hz. Valid input are 50 and 60
%
% Returns
% -------
%     harmonic_subgroups: structure with handles harmonic_subgroups.amplitude and
%           harmonic_subgroups.harmonic
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

sg_pd = py.mhkit.power.quality.harmonic_subgroups(data_pd,grid_freq);

vals=double(py.array.array('d',py.numpy.nditer(sg_pd.values)));
sha=cell(sg_pd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[y,x]).';


harmonic_subgroups.amplitude=vals;
harmonic_subgroups.harmonic = double(py.array.array('d',py.numpy.nditer(sg_pd.index)));
