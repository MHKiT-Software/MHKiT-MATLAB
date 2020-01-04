function J=energy_flux(S,h,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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
%    h: float
%         Water depth (m)
%
%    rho: float (optional)
%         water density (kg/m^3)
%    g: float (optional)
%         gravitational acceleration (m/s^2)
%
%     NOTE: In matlab, if you set one optional parameter, you must set
%     both, rho first, then g
%         
% Returns
% -------
%     J: double
%         Omni-directional wave energy flux (W/m)
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
            
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(uint32(S.frequency(:,1)),li,int32(x(2)));
        elseif x(2)==1
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(uint32(S.frequency),py.numpy.array(S.spectrum),int32(x(2)));
        end
    else
        ME = MException('MATLAB:energy_flux','S needs to be structure or a Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end


if nargin == 4 
    J=py.mhkit.wave.resource.energy_flux(S,h,pyargs('rho',varargin{1},'g',varargin{2}));
elseif nargin == 2
    J=py.mhkit.wave.resource.energy_flux(S,h);
else
    ME = MException('MATLAB:energy_flux','incorrect number of arguments');
        throw(ME);
end


J=double(J.values);
