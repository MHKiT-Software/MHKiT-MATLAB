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



py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

if (isa(S,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(S)==1)
        x=size(S.spectrum);
        li=py.list();
        if x(2)>1 
            for i = 1:x(2)
                app=py.list(S.spectrum(:,i));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
            end
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency(:,1),li,x(2));
        elseif x(2)==1
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency,py.numpy.array(S.spectrum),x(2));
        end
    else
        ME = MException('MATLAB:frequency_moment','S needs to be a Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

if nargin == 3
    m=py.mhkit.wave.resource.frequency_moment(S,int32(N),pyargs('frequency_bins',py.numpy.array(varargin{1})));
elseif nargin == 2
    m=py.mhkit.wave.resource.frequency_moment(S,int32(N));
else
    ME = MException('MATLAB:frequency_moment','Incorrect number of arguments');
        throw(ME);
end

m=double(m.values);