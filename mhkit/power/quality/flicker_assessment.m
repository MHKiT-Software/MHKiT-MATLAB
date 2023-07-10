function [alpha0,freq,alpha_m,u0,u_fic] = flicker_assessment(Sr,Un,SCR,fg, ...
        u_m,i_m,method,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduct flicker assessment according to IECTS-30
%   Workflow of MATLAB functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Construct the fictitious grid:
[Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
%% 2. Calculate ideal voltage (u0) from measured voltage (u_m)
% 2.1 Calculate the fundamental frequency and alpha_0 of u_m:
% STFT or ZCD
[alpha0,freq] = calc_alpha0_freq(u_m,method,opts);
% 2.2 Calculate the electrical angle (alpha_m):
alpha_m = calc_electrical_angle(freq,alpha0);
% 2.3 Calculate ideal voltage (u0):
u0 = calc_ideal_voltage(Un,alpha_m);
%% 3. Calculate simulated voltage (u_fic)
u_fic = calc_simulated_voltage(u0,i_m,Rfic,Lfic);
end

