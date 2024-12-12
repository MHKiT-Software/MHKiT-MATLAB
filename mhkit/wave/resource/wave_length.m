function l=wave_length(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave length from wave number
% To compute: 2*pi/wavenumber
%
% Parameters
% ------------
%    k: wave number (1/m)
%       intiger, double, or vector
%
% Returns
% -------
%    l: double or array
%         Wave length [m] indexed by frequency
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=double(py.array.array('d',py.numpy.nditer(py.mhkit.wave.resource.wave_length(py.numpy.array(k)))));

