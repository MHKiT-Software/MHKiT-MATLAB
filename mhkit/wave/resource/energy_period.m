function Te=energy_period(S,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Parameters
% ------------
%   S: Spectral Density (m^2/Hz)
%       Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,spectra)
%
%       OR
%
%       structure of form:
%           S.spectrum: Spectral Density (m^2/Hz)
%
%           S.type: String of the spectra type, i.e. Bretschneider,
%           time series, date stamp etc.
%
%           S.frequency: frequency (Hz)
%
%     frequency_bins: vector (optional)
%       Bin widths for frequency of S. Required for unevenly sized bins
%
%
% Returns
% ---------
%    Te: float
%        Wave energy Period (s)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign feq_bin
if nargin == 2
    freq_bins = py.numpy.array(varargin{1});
elseif nargin == 1
    freq_bins = py.None;
else
    ME = MException('MATLAB:energy_period','Incorrect number of input arguments');
        throw(ME);
end

S_py = typecast_spectra_to_mhkit_python(S);

Te = py.mhkit.wave.resource.energy_period(S_py, pyargs('frequency_bins',freq_bins));

Te = typecast_from_mhkit_python(Te).data;
