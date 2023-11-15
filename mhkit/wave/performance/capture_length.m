function L=capture_length(Power,J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the capture length (often called capture width).
%
% Parameters
% ------------
%     P: array or vector
%         Power [W]
%
%     J: array or vector
%         Omnidirectional wave energy flux [W/m]
%
% Returns
% ---------
%     L: vector
%         Capture length [m]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');

Power=py.numpy.array(Power);
J=py.numpy.array(J);


L=double(py.mhkit.wave.performance.capture_length(Power,J)).';

