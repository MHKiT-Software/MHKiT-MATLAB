function [Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates electrical angle of the fundamental of the measured
%   voltage. According to IEC Standard 62600-30 formula(3).
%
%   
% Parameters
% -----------
%   Sr: double
%       rated apparent power of the marine energy converter (VA). 
%   SCR: double
%       short-circuit ratio, SCR = S_kfic/Sr
%       To obtain simulated voltage fluctuations within the flicker meter 
%       range, an SCR of [20,50] is recommended. (ref: IECTS62600-30(ed1.0)
%       pg 24)
%   Un: double
%       RMS value of the nominal voltage of the grid 
%   fg: double
%       nominal grid frequency, 50Hz or 60Hz
% 
% Returns
% ------- 
%   Rfic: double array of size (4) 
%       fictitious grid resistance (Ohm) for the impedance phase angle (Phi_k)
%       equals to 30, 50, 70, and 85 degree.
%   Lfic: double array of size (4)
%       fictitious grid inductance (H) for Phi_k equals to 30, 50, 70, 
%       and 85 degree.
% Note
% -------
%   Phi_k: network impedance phase angle (degree)
%       The flicker index should be reported for Phi_k=30, 50, 70, 85
%       (degree). (refs: IECTS62600-30(ed1.0) Annex A Table A.6)
%   Formulas from IECTS62600-30
%       SCR,short-circuit ratio:
%           SCR = S_kfic/Sr (pg24)
%       S_kfic, short-circuit apparent power (VA):
%           S_kfic = Un^2/sqrt(Rfic^2+X_fic^2) (Eq5) 
%       Phi_k, network impedance phase angle:
%           tan(Phi_k) = X_fic/Rfic = 2*pi*fg*Lfic/Rfic (Eq4)
%       X_fic, reactance of fictitious grid (Ohm)
%   Therefore, we have:
%       sqrt{Rfic^2+[Rfic*tan(Phi_k)]^2} = Un^2/S_k,fic -->>
%
%       Rfic = Un^2/{S_kfic*sqrt[1+tan(Phi_k)^2]} 
%   and
%       Lfic = Rfic*tan(Phi_k)/(2pi*fg)
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % impedance angles in degrees
    Phi_k = [30,50,70,85]; 
    S_kfic = SCR*Sr;
    Rfic = Un.^2./(S_kfic.*sqrt(1+(tand(Phi_k)).^2));
    X_fic = Rfic.*tand(Phi_k);
    Lfic = X_fic/(2*pi*fg);
end