function k=wave_number(f,h,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave number
%
% Parameters
% ------------
%    f: frequency (Hz)
%           vector or numpy array
%    h: float
%         Water depth (m)
%
%    rho: float (optional)
%         water density (kg/m^3)
%         to call: wave_number(f,h,"rho",rho)
%
%    g: float (optional)
%         gravitational acceleration (m/s^2)
%         to call: wave_number(f,h,"g",g)
%
%
%
%
% Returns
% -------
%     k: structure
%
%
%         k.values: wave number
%
%         k.frequency: frequency [Hz]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    f
    h
    options.rho = 1025;
    options.g = 9.80665;
end

f=py.numpy.array(f);

kdf=py.mhkit.wave.resource.wave_number(f,h,pyargs('rho',options.rho,'g',options.rho));

k.values=double(kdf.values);
k.frequency=double(kd.index);
