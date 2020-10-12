function J=energy_flux(S,h,options)

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
%         to call: energy_flux(S,h,"rho",rho)
%
%    g: float (optional)
%         gravitational acceleration (m/s^2)
%         to call: energy_flux(S,h,"g",g)
%
% Returns
% -------
%     J: double
%         Omni-directional wave energy flux (W/m)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S 
    h
    options.rho = 1025;
    options.g = 9.80665;
    
end

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
            
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency(:,1),li,int32(x(2)));
        elseif x(2)==1
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency,py.numpy.array(S.spectrum),int32(x(2)));
        end
    else
        ME = MException('MATLAB:energy_flux','S needs to be structure or a Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end



J=py.mhkit.wave.resource.energy_flux(S,h,pyargs('rho',options.rho,'g',options.g));



J=double(J.values);
