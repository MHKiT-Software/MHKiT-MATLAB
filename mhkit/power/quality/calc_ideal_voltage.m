function u0=calc_ideal_voltage(Un,alpha_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the ideal phase-to-neutral voltage source (V). According to 
% IEC standard 62600-30(ed1.0) formula (2).
%
%   The ideal voltage should: 1) be without any fluctuations; 2) have the 
% same electrical angle (alpha_m) as the fundamental of the measured
% voltage (u_m).
%   
% Parameters
% -----------
%   Un: double 
%       RMS value of the nomianal voltage of the grid (V).
%   alpha_m: double array
%       the electrical angle (radian) of the fundamental of the measured 
%       voltage. 
% 
% Returns
% -------
%   u0: struct()
%       .time: time for each time step; 
%       .data: ideal phase-to-neutral voltage source with the instantaneous
%       value u0(t).
% Note
% -------
% 1. According to the IEC standard 62600-30(ed1.0) formula (2):
%   u0(t) = sqrt(2/3)*Un*sin(alpha_m(t)).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IEC standard 62600-30(ed1.0) formula (2)
    u0 = sqrt(2/3)*Un*sin(alpha_m);
end