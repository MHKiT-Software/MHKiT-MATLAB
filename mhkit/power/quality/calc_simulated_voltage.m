function u_fic=calc_simulated_voltage(u0,i_m,Rfic,Lfic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the simulated voltage at a fictitious grid according to 
%   IEC standard 62600-30(ed1.0) Eq (1).
%   
% Parameters
% -----------
%   u0: double array of size (ntime)
%       Ideal phase-to-neutral voltage source u0(t) (V).
%   i_m: struct()
%       i_m(t) is the measured instantaneous current (A)
%       .time: time at each measurement
%       .data: array of size (ntime, 4) measured instantaneous current (A)
%       at impedance phase angles Phi_k = 30, 50, 70, 85.
%   Rfic: double array of size (4) 
%       Fictitious grid resistance (Ohm) for Phi_k = 30, 50, 70, 85.
%   Lfic: double array of size (4)
%       Fictitious grid inductance (H) for Phi_k = 30, 50, 70, 85.
% 
% Returns
% -------
%   u_fic: double array (ntime,4) 
%       Instantaneous phase-to-neutral voltage simulated at fictitious 
%       grid (V) for Phi_k = 30, 50, 70, 85
% Note
% ------- 
%   1. ufic(t) = u0(t)+Rfic*im(t)+Lfic*dim(t)/dt (Eq1, IECTS62600-30)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % 1. calculate dim(t)/dt using diff()
    delta_im = reshape(diff(i_m.data),[],4);
    delta_t = reshape(diff(i_m.time),[],1);
    dim_dt = delta_im./delta_t;
    dim_dt(end+1,:) = dim_dt(end,:);
    % Can also consider using: Center-Euler, Euler-Forward, OR 
    % Euler-Backward to calculate dim/dt
    % 2. calculate u_fic(t) according to Eq (1) of IECTS62600-30
    u_fic = u0 + i_m.data.*Rfic+dim_dt.*Lfic;

end