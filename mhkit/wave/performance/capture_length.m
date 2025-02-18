function L = capture_length(Power, J)

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



Power = py.numpy.array(Power);
J = py.numpy.array(J);

L = py.mhkit.wave.performance.capture_length(Power, J);

L = typecast_from_mhkit_python(L).data;
