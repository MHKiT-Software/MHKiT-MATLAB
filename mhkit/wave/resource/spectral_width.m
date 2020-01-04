function v=spectral_width(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave spectral width from spectra
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
% Returns
% ---------
%    e0: float
%        Spectral Width
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
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(S.frequency(:,1)),li,int32(x(2)));
        elseif x(2)==1
            S=py.mhkit+python_utils.pandas_dataframe.spectra_to_pandas(double(S.frequency),py.numpy.array(S.spectrum),int32(x(2)));
        end
       
    else
        ME = MException('MATLAB:significant_wave_height','S needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

e0=py.mhkit.wave.resource.spectral_width(S);
v=double(e0.values);