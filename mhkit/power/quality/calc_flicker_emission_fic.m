function P_stfic = calc_flicker_emission_fic(u_fic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates the flicker emission value P_st,fic on the fictitious grid
%   for each 10 min time-series, which is then used to calculate the 
%   flicker coefficient (c). The calculation is according to the 
%   calculation of the short-term flicker severity (Pst) in the IEC 
%   61000-4-15(ed2.0) ***.
%
% Parameters
% -----------
%   u_fic: struct()
%       .time: time for each time step; 
%       .data: voltage time series of the u_fic(t), instantaneous 
%       phase-to-neutral voltage (V) simulated on the fictitious grid.
% 
% Returns
% -------
%   P_stfic: double array of size dtime/10min 
%       One flicker emission value on the fictitious grid for each 10 min 
%       time-series of u_fic(t).
%           
% Note
% -------
% 1. Unless otherwise specified, the Pst evaluation time is 10min. For the
%   purpose of power quality surveys and studies, other time intervals may 
%   be used, and should be defined in the index. For example a 1 min 
%   interval should be written as P_st,1m.
% 2. Tshort: short term interval for the Pst evaluation, unless otherwise
%   specified, Tshort is 10 min.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate Pst according to IEC61600-4-15
    P_stfic = struct();
    
end