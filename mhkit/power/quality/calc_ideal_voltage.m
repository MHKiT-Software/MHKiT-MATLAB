function u0=calc_ideal_voltage(Un,alpha_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the ideal phase-to-neutral voltage source (V) u0(t)
% according to IECTS 62600-30(ed1.0) Eq (2).
%
%   u0(t) should:
%   1) be without any fluctuations;
%   2) have the same electrical angle (alpha_m) as the fundamental of the
%   measured voltage (u_m).
%
% Parameters
% -----------
%   Un: double
%       RMS value of the nominal voltage of the grid (V).
%   alpha_m: double array (ntime)
%       Electrical angle of the fundamental component of u_m(t)
% Returns
% -------
%   u0: double array (ntime)
%       Ideal phase-to-neutral voltage source (V).
% Note
% -------
% 1. According to the IECTS-62600-30(ed1.0) Eq (2):
%   u0(t) = sqrt(2/3)*Un*sin(alpha_m(t)).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IECTS 62600-30(ed1.0) Eq (2)
    u0 = sqrt(2/3)*Un*sin(alpha_m);
end

