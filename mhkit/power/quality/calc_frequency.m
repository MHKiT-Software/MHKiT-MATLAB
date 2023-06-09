function freq = calc_frequency(u_m,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the frequency of the measured voltage u_m(t).
%   
% Parameters
% -----------
%   u_m: struct()
%       .time: time for each time step; 
%       .data: measured voltage with instantaneous value u_m(t) (V).
%   method: ???
%       option of methods to use to calculate alpha_m, DFT or ZCD.
% 
% Returns
% -------
%   freq: struct()
%       .time: time for each time step; 
%       .data: freq(t), the frequency for measured voltage.
% Note
% -------
%   Two methods can be used to calculate the frequency:
%   1) STFT/DFT: Short Time Fourier Transform
%   2) ZCD: Zero Crossing Detection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate frequency using DFT or ZCD
    freq = struct();
end