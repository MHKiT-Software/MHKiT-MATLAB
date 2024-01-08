function J=energy_flux(S,h,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Parameters
% ------------
%    S: Spectral Density (m^2/Hz)
%       structure of form:
%           S.spectrum: Spectral Density (m^2/Hz)
%
%           S.type: String of the spectra type, i.e. Bretschneider,
%           time series, date stamp etc.
%
%           S.frequency: frequency (Hz)
%
%    h: double
%         Water depth (m)
%
%    deep: logical (optional)
%         If True use the deep water approximation. Default False. When
%         False a depth check is run to check for shallow water. The ratio
%         of the shallow water regime can be changed using the ratio
%         keyword.
%         to call: energy_flux(S,h,"deep",py.True)
%
%    rho: double (optional)
%         water density (kg/m^3)
%         to call: energy_flux(S,h,"rho",rho)
%
%    g: double (optional)
%         gravitational acceleration (m/s^2)
%         to call: energy_flux(S,h,"g",g)
%
%    ratio: double or int (optional)
%         Only applied if depth=False. If h/l > ratio,
%         water depth will be set to deep. Default ratio = 2.
%         to call: energy_flux(S,h,"ratio",1.5)
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
    options.deep = py.False;
    options.rho = 1025;
    options.g = 9.80665;
    options.ratio = 2;

end

py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

% convert to dataframe
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


% calculate energy flux
J=py.mhkit.wave.resource.energy_flux(S,h,pyargs('rho',options.rho,'g',options.g,...
    'deep',options.deep, 'ratio',options.ratio));



J=double(J.values);

