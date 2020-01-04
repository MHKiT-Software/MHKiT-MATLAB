function clm=capture_length_matrix(Hm0,Te,L,statistic,Hm0_bins,Te_bins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Generates a capture length matrix for a given statistic
%    
%     Note that IEC/TS 62600-100 requires capture length matrices for 
%     the mean, std, count, min, and max.
%     
% Parameters
% ------------
%     Hm0: numpy array or vector
%         Significant wave height from spectra [m]
%
%     Te: numpy array or vector
%         Energy period from spectra [s]
%
%     L : numpy array or vector
%         Capture length [m]
%
%     statistic: string
%         Statistic for each bin, options include: 'mean', 'std', 'median', 
%         'count', 'sum', 'min', 'max', and 'frequency'.  Note that 'std' uses 
%         a degree of freedom of 1 in accordance with IEC/TS 62600-100. or
%         a callable function of python type
%
%     Hm0_bins: numpy array or vector
%         Bin centers for Hm0 [m]
%
%     Te_bins: numpy array or vector
%         Bin centers for Te [s]
%         
% Returns
% ---------
%     clm: structure
%
%
%        clm.values
%
%        clm.stat
%
%        clm.Hm0_bins
%
%        clm.Te_bins
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');

Hm0=py.numpy.array(Hm0);
Te=py.numpy.array(Te);
L=py.numpy.array(L);
Hm0_bins=py.numpy.array(Hm0_bins);
Te_bins=py.numpy.array(Te_bins);

LM=py.mhkit.wave.performance.capture_length_matrix(Hm0,Te,L,statistic,Hm0_bins,Te_bins);
vals=double(py.array.array('d',py.numpy.nditer(LM.values)));
sha=cell(LM.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[y,x]);
vals=transpose(vals);

clm.values=vals;
clm.stat=statistic;
clm.Hm0_bins=double(Hm0_bins);
clm.Te_bins=double(Te_bins);
 
