function WEFM=wave_energy_flux_matrix(Hm0,Te,J,statistic,Hm0_bins,Te_bins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Generates a wave eneergy flux matrix for a given statistic
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
%     J : numpy array or vector
%         wave energy flux from spectra [W/m]
%
%     statistic: string
%         Statistic for each bin, options include: 'mean', 'std', 'median',
%         'count', 'sum', 'min', 'max', and 'frequency'.  Note that 'std' uses
%         a degree of freedom of 1 in accordance with IEC/TS 62600-100.
%
%     Hm0_bins: numpy array or vector
%         Bin centers for Hm0 [m]
%
%     Te_bins: numpy array or vector
%         Bin centers for Te [s]
%
% Returns
% ---------
%     WEFM: Structure
%
%
%       WEFM.values
%
%       WEFM.stat
%
%       WEFM.Hm0_bins
%
%       WEFM.Te_bins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');

Hm0=py.numpy.array(Hm0);
Te=py.numpy.array(Te);
J=py.numpy.array(J);
Hm0_bins=py.numpy.array(Hm0_bins);
Te_bins=py.numpy.array(Te_bins);

JM=py.mhkit.wave.performance.wave_energy_flux_matrix(Hm0,Te,J,statistic,Hm0_bins,Te_bins);
vals=double(py.array.array('d',py.numpy.nditer(JM.values)));
sha=cell(JM.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[y,x]);
vals=transpose(vals);

WEFM.values=vals;
WEFM.stat=statistic;
WEFM.Hm0_bins=double(Hm0_bins);
WEFM.Te_bins=double(Te_bins);

