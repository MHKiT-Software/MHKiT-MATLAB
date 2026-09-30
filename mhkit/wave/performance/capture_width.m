function CW = capture_width(P, J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the capture width (sometimes called capture length)
%
% CW = P / J
%
% Parameters
% ------------
% P : numeric [W]
%   Power (scalar, vector, or array)
% J : numeric [W/m]
%   Omnidirectional wave energy flux (same shape as P)
%
% Returns
% ---------
% CW : numeric [m]
%   Capture width (same shape as P and J)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    P {mustBeNumeric}
    J {mustBeNumeric}
end

arguments (Output)
    CW {mustBeNumeric}
end

CW = P ./ J;

end
