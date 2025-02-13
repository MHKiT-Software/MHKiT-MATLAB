function m=frequency_moment(S,N,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the Nth frequency moment of the spectrum
%
% Parameters
% ------------
%    S: Spectral Density (m^2/Hz)
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
%    N: int
%       Moment (0 for 0th, 1 for 1st ....)
%
%    frequency_bins: vector (optional)
%       Bin widths for frequency of S. Required for unevenly sized bins
%
% Returns
% ---------
%    m: double
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_py = typecast_spectra_to_mhkit_python(S);

if nargin == 3
    m = py.mhkit.wave.resource.frequency_moment(S_py, int32(N),pyargs('frequency_bins',py.numpy.array(varargin{1})));
elseif nargin == 2
    m = py.mhkit.wave.resource.frequency_moment(S_py, int32(N));
else
    ME = MException('MATLAB:frequency_moment','Incorrect number of arguments');
        throw(ME);
end

m = typecast_from_mhkit_python(m).data;
