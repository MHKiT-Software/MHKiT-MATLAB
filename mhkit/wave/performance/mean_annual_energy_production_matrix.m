function maep=mean_annual_energy_production_matrix(LM,JM,frequency)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculates mean annual energy production (MAEP) from matrix data
%     along with data frequency in each bin
%
% Parameters
% ------------
%     LM: Capture Length
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,L)
%
%        OR
%
%        structure of form:
%
%           LM.values
%
%           LM.stat
%
%           LM.Hm0_bins
%
%           LM.Te_bins
%
%
%     JM: Wave Energy Flux
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,J)
%
%        OR
%
%        structure of form:
%
%           JM.values: Wave energy flux matrix
%
%           JM.Hm0_bins
%
%           JM.Te_bins
%
%     frequency: Data frequency for each bin.
%         Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,frequency)
%
%        OR
%
%        structure of form:
%
%           frequency.values
%
%           frequency.Hm0_bins
%
%           frequency.Te_bins
%
% Returns
% ---------
%     maep: float
%         Mean annual energy production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LM_py = py.pandas.DataFrame(py.numpy.array(LM.values));
JM_py = py.pandas.DataFrame(py.numpy.array(JM.values));
freq_py = py.pandas.DataFrame(py.numpy.array(frequency.values));

maep = py.mhkit.wave.performance.mean_annual_energy_production_matrix(LM_py,JM_py,freq_py);
maep = double(maep);
