function u0=calc_ideal_voltage(Un,u_m,alpha_0,method)
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
%   u_m: struct()
%       .time: time for each time step; 
%       .data: measured voltage with instantaneous value u_m(t) (V).
%   alpha_0: double 
%       the electrical angle at t=0 (radian). 
%   method: ???
%       option of methods to use to calculate alpha_m, DFT or ZCD.
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
%   alpha_m(t) = 2pi*integral(freq(t),[0,t])+alpha_0
% 3. freq(t): frequency of the measured u_m (that may vary over time).
%   freq = calc_frequency(u_m,method);
%   Two methods can be used to calculate the frequency:
%   1) STFT/DFT: Short Time Fourier Transform
%   2) ZCD: Zero Crossing Detection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate frequency using u_m and corresponding method:
    freq = calc_frequency(u_m,method);

    % calculate alpha_m:
    % use trapz() to calculate integral of freq(t) from [0,t]
    alpha_m = 2*pi*trapz(freq,u_m.time)+alpha_0;

    %IEC standard 62600-30(ed1.0) formula (2)
    u0 = sqrt(2/3)*Un*sin(alpha_m);
end