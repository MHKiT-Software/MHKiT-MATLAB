function c=calc_flicker_coefficient(P_stfic,S_kfic,Sr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the flicker coefficient c(Phi_k) for each set of the 10 min 
%   measured voltage and current time-series, using the calculated 
%   flicker emission values by applying Eq 6 from IECTS-62600-30. 
%   
%   
% Parameters
% -----------
%   P_stfic: double 
%       short-term flicker severity (Pst) (flicker emission value) for the
%       fictitious grid. 
%   S_kfic: double
%       short-circuit apparent power (VA) of the fictitious grid.
%       S_kfic = Un^2/sqrt(Rfic^2+X_fic^2) (Eq5) 
%   SCR,short-circuit ratio:
%           SCR = S_kfic/Sr (pg24)
%   Sr: double
%       rated apparent power (VA) of the marine energy converter unit. 
% 
% Returns
% -------
%   c: double array of size (4) 
%       flicker coefficient for the impedance phase angle (Phi_k)
%       equals to 30, 50, 70, and 85 degree.
% Note
% -------
% 1. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate flicker coefficient:
    c = P_stfic*S_kfic/Sr;
end