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



% Flatten arrays to 1D to handle both row and column vectors
Power = py.numpy.array(Power).flatten();
J = py.numpy.array(J).flatten();

L = py.mhkit.wave.performance.capture_width(Power, J);

L = typecast_from_mhkit_python(L).data;
