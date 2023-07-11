function [alpha0,freq,alpha_m,u0,u_fic] = flicker_ufic(Sr,Un,SCR,fg, ...
        u_m,i_m,method,methodopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduct flicker assessment according to IECTS62600-30 and
%   IECTS61400-21-1.
%   A workflow of MATLAB functions up to calculation of u_fic. 
% 
% Parameters
% -----------
%   Sr: double 
%       Rated apparent power of the marine energy converter (VA).
%   Un: double 
%       RMS value of the nomianal voltage of the grid (V).
%   SCR: double 
%       Short-circuit ratio, SCR = S_kfic/Sr
%   fg: double 
%       Nominal grid frequency, 60Hz or 50Hz
%   u_m: struct() 
%       u_m(t) is the measured instantaneous voltage (V)
%       .time: time at each measurement (s)
%       .data: array of size (ntime) measured instantaneous voltage (V).
%   i_m: struct()
%       i_m(t) is the measured instantaneous current (A)
%       .time: time at each measurement
%       .data: array of size (ntime, 4) measured instantaneous current (A)
%       at impedance phase angles = 30, 50, 70, 85.
%   method: string 
%       Method used to calculate fundamental frequency of u_m, 'zcd' or 
%       'stft'
%   methodopts: cell array 
%       Name-value arguments for STFT method, if using 'ZCD' method, can 
%       set to {}.    
% 
% Returns
% -------
%   alpha0: double 
%       Electrical angle (radians) at t=0.
%   freq: struct()
%       .time: time at each measurement time step (s) 
%       .data: double array of size (ntime), freq(t), the fundamental 
%       frequency calculated for u_m(t).
%   alpha_m: double array (ntime)
%       Electrical angle (alpha_m(t)) of the fundamental component of u_m(t)
%   u0: double array (ntime)
%       Ideal phase-to-neutral voltage source u0(t) (V)
%   u_fic: double array (ntime,4) 
%       Instantaneous phase-to-neutral voltage simulated at fictitious 
%       grid (V) for impedence angles = 30, 50, 70, 85
% Note
% -------
% Step 1. Construct the fictitious grid by calculating Rfic and Lfic.
%       [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg)
%
% Step 2. Calculate ideal voltahe source (u0) from measured voltage (u_m).
%   Step 2.1 Calculate the fundamental frequency and alpha_0 of u_m 
%       [alpha0,freq] = calc_alpha0_freq(u_m,method,methodopts) 
%   Step 2.2 Calculate the electrical angle (alpha_m) 
%       alpha_m = calc_electrical_angle(freq,alpha0)
%   Step 2.3 Calculate the ideal voltage (u0) 
%       u0=calc_ideal_voltage(Un,alpha_m)
% 
% Step 3. Calculate simulated voltage (u_fic) for the fictitious grid 
%       u_fic=calc_simulated_voltage(u0,i_m,Rfic,Lfic)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Construct the fictitious grid:
[Rfic,Lfic] = calc_Rfic_Lfic(Sr, SCR, Un, fg);
%% Step 2. Calculate ideal voltage (u0) from measured voltage (u_m)
% Step 2.1 Calculate the fundamental frequency and alpha_0:
% opt1: method ZCD: method = 'ZCD'; methodopts={};
% opt2: method STFT: method = 'stft'; methodopts = {'Window',rectwin(M),...
% 'OverlapLength',L,'FFTLength',128,'FrequencyRange','onesided'};
[alpha0,freq] = calc_alpha0_freq(u_m,method,methodopts);
% Step 2.2 Calculate the electrical angle (alpha_m):
alpha_m = calc_electrical_angle(freq,alpha0);
% Step 2.3 Calculate ideal voltage (u0):
u0 = calc_ideal_voltage(Un,alpha_m);
%% Step 3. Calculate simulated voltage (u_fic)
u_fic = calc_simulated_voltage(u0,i_m,Rfic,Lfic);
end

