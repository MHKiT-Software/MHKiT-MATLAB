function u0=calc_ideal_voltage(Un,alpha0,freq)
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
%   alpha0: double
%       the electrical angle at t=0 (radian). 
%   freq: struct()
%       .time: time instants corresponding to the centers of the data 
%              segments used to capture the power spectrum estimates. 
%       .data: freq(t), the fundamental frequency calculated from the 
%              STFT for u_m. Frequency of the measured u_m (that may vary 
%              over time).
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
%   u0(t) = sqrt(2/3)*Un*sin(alpha_m(t)), where:
%
% 2. alpha_m(t): electrical angle (radian) of the fundamental of the 
%   measured voltage, according to the IEC Standard 62600-30 formula(3):
%   alpha_m(t) = 2pi*integral(freq(t)dt)+alpha0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate alpha_m:
    alpha_m = 2*pi*cumtrapz(freq.time,freq.data)+alpha0;
    %IEC standard 62600-30(ed1.0) formula (2)
    u0 = sqrt(2/3)*Un*sin(alpha_m);
end