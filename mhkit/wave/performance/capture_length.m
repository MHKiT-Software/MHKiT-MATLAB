function L = capture_length(P, J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Deprecated alias for capture_width.
%
% IEC TS 62600-100 Ed. 2.0 replaces "capture length" with "capture
% width". This function will be removed in MHKiT-MATLAB v1.3; use
% capture_width instead.
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
% L : numeric [m]
%   Capture length (same shape as P and J)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    P {mustBeNumeric}
    J {mustBeNumeric}
end

arguments (Output)
    L {mustBeNumeric}
end

warning('MHKiT:capture_length:DeprecatedFunction', ...
    ['IEC TS 62600-100 Ed. 2.0 replaces "capture length" with "capture ' ...
    'width". capture_length will be removed in MHKiT-MATLAB v1.3. ' ...
    'Use capture_width instead.']);

L = capture_width(P, J);

end
