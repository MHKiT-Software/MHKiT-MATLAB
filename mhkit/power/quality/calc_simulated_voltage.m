function u_fic=calc_simulated_voltage(u0,i_m,Rfic,Lfic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the simulated voltage at a fictitious grid according to 
%   IEC standard 62600-30(ed1.0) formula (1).
%   
% Parameters
% -----------
%   u0: struct()
%       .time: time for each time step (=im.time).
%       .data: ideal phase-to-neutral voltage source with the instantaneous 
%       value u0(t) (V).
%   i_m: struct()
%       .time: time at each measurement
%       .data: array of size (ntime) measured instantaneous current (A).
%   Rfic: double array of size (4) 
%       fictitious grid resistance for the impedance phase angle (Phi_k)
%       equals to 30, 50, 70, and 85 degree (Ohm).
%   Lfic: double array of size (4)
%       fictitious grid inductance for the impedance phase angle (Phi_k)
%       equals to 30, 50, 70, and 85 degree (H).
% 
% Returns
% -------
%   u_fic: instantaneous phase-to-neutral voltage simulated at 
%       a fictitious grid (V).
% Note
% ------- 
%   1. ufic(t) = u0(t)+Rfic*im(t)+Lfic*dim(t)/dt (Eq1, IECTS62600-30)
%   2. Prious steps: calc_fictitious_grid & calc_ideal_voltage
%       [Rfic, Lfic] = calc_fictitious_grid(Sr, SCR, Un, fg) 
%       u0 = calc_ideal_voltage(Un, u_m, alpha_0, method)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    % calculate dim(t)/dt
    % use diff() function in MATLAB
    delta_im = reshape(diff(i_m.data),[],1);
    delta_t = reshape(diff(i_m.time),[],1);
    dim_dt = delta_im./delta_t;
    dim_dt(end+1) = dim_dt(end);
    % Can also consider using: 
    % Center-Euler, Euler-Forward, OR Euler-Backward
    % to calculate dim/dt
    % calculate u_fic(t)
    %u_fic = reshape(u0.data,[],1) + Rfic.*i_m.data+Lfic.*dim_dt;
    disp(size(u0.data));disp(size(Rfic))
    u_fic = u0.data + Rfic'.*i_m.data+Lfic'.*dim_dt';
    


end