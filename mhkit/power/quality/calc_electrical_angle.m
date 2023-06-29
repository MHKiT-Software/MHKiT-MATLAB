function alpha_m = calc_electrical_angle(freq,alpha0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the electrical angle alpha_m(t) of the fundamental of the 
% measured voltage. According to IEC standard 62600-30(ed1.0) formula (3).
%
%   
% Parameters
% -----------
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
%   alpha_m: double array
%       alpha_m(t) is the electrical angle (radian) for the fundamental of 
%       the measured voltage (u_m) at each time step.
% Note
% -------
% 1. According to the IEC Standard 62600-30 formula(3):
%   alpha_m(t) = 2pi*integral(freq(t)dt)+alpha0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculate alpha_m:
    alpha_m = 2*pi*cumtrapz(freq.time,freq.data)+alpha0;
end