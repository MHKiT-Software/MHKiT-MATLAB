function coef_flicker=calc_flicker_coefficient(P_stfic,S_kfic,Sr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the flicker coefficient c(Phi_k) for each set of the 10 min
%   measured voltage and current time-series, using the calculated
%   flicker emission values using IECTS-62600-30 Eq (6).
%
%
% Parameters
% -----------
%   P_stfic: double
%       short-term flicker severity (Pst) (flicker emission value) on the
%       fictitious grid.
%   S_kfic: double
%       Short-circuit apparent power (VA) of the fictitious grid.
%       S_kfic = Un^2/sqrt(Rfic^2+X_fic^2) (Eq5)
%   SCR: double
%       Short-circuit ratio:
%           SCR = S_kfic/Sr (pg24)
%   Sr: double
%       Rated apparent power (VA) of the marine energy converter unit.
%
% Returns
% -------
%   coef_flicker: double
%       Flicker coefficient for continuous operation
% Note
% -------
% 1. According to IECTS62600-30 Eq (6):
%       c(Phi_k) = P_stfic*S_kfic/Sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % IECTS62600-30 Eq (6):
    coef_flicker = P_stfic*S_kfic/Sr;
end

