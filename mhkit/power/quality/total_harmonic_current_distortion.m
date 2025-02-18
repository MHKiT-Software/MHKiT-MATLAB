function THCD=total_harmonic_current_distortion(harmonic_subgroups,rated_current)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the total harmonic current distortion (THC) based on IEC/TS 62600-30
%
% Parameters
% -----------
%     harmonic_subgroups: structure with handles- harmonic_subgroups.amplitude and harmonic_subgroups.harmonic
%         Subgrouped current harmonics amplitude indexed by harmonic order
%
%     rated_current: double
%         Rated current of the energy device in Amps
%
% Returns
% -------
%     THCD: double
%         Total harmonic current distortion
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = harmonic_subgroups.amplitude;

dsize=size(data);
li=py.list();
if dsize(2)>1
   for i = 1:dsize(2)
      app=py.list(data(:,i));
      li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

   end
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(harmonic_subgroups.harmonic,li,int32(dsize(2)));
elseif dsize(2)==1
   data_pd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(harmonic_subgroups.harmonic,py.numpy.array(data),int32(dsize(2)));
end

thcd_pd = py.mhkit.power.quality.total_harmonic_current_distortion(data_pd);
THCD = double(py.array.array('d',py.numpy.nditer(thcd_pd.values)));
