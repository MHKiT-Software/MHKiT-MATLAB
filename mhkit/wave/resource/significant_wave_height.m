function Hm0=significant_wave_height(S,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave height from spectra
%
% Parameters
% ------------
%     S: Spectral Density (m^2/Hz)
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
% Returns
% ---------
%     Hm0: double 
%         Significant Wave Height (m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

if nargin == 2
    freq_bins = py.numpy.array(varargin{1});
elseif nargin == 1
    freq_bins = py.None;
else
    ME = MException('MATLAB:significant_wave_height','Incorrect number of input arguments');
        throw(ME);
end
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
        ME = MException('MATLAB:significant_wave_height','S needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

Hm0=py.mhkit.wave.resource.significant_wave_height(S,pyargs('frequency_bins',freq_bins));
Hm0=double(Hm0.values);
