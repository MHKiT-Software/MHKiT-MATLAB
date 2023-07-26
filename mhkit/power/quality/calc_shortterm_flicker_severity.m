function Pst = calc_shortterm_flicker_severity(P)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates short-term flicker severity P_st according to the IEC 
%   61000-4-15(ed2.0) section 5.7.2 Short-term flicker evaluation.
%
% Parameters
% -----------
%   P: struct()
%       contains fields of flicker levels exceeded for *.* percent of the
%       time. 
%       P.p0p7, P.p1, P.p1p5: flicker levels exceeded for 0.7%, 1%, and 
%       1.5% of the time. 
%       P.p2p2, P.p3, P.p4: flicker levels exceeded for 2.2%, 3%, and 4%
%       of the time.
%       P.p6, P.p8, P.p10, P.p13, P.p17: flicker levels exceeded for 6%,
%       8%, 10%, 13%, 17% of the time.
%       P.p30, P.p50, P.p80: flicker levels exceeded for 30%, 50%, and 80%
%       of the time.
%       
% 
% Returns
% -------
%   P_st: double array of size dtime/10min 
%       One flicker emission value on the fictitious grid for each 10 min 
%       time-series of u_fic(t).
%           
% Note
% -------
% 1. Unless otherwise specified, the Pst evaluation time is 10min. For the
%   purpose of power quality surveys and studies, other time intervals may 
%   be used, and should be defined in the index. For example a 1 min 
%   interval should be written as P_st,1m.
% 2. The suffix '_s' indicates the smoothed values and obtained following
% the equations specificed in section 5.7.2:
%   P_50s = (P_30+P_50+P_80)/3.;
%   P_10s = (P_6+P_8+P_10+P_13+P_17)/5.;
%   P_3s = (P_2.2+P_3+P_4)/3.;
%   P_1s = (P_0.7+P_1+P_1.5)/3.;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check input:
    names_need = {'p30','p50','p80','p6','p8','p10','p13','p17',...
        'p2p2','p3','p4','p0p7','p1','p1p5','p0p1'};
    for i=1:length(names_need)
        if ~isfield(P,names_need{i})
            ME = MException('MATLAB:calc_shortterm_flicker_severity',...
            'invalid handles in structure, must contain x.data & x.time');
            throw(ME);
        end
    end
    P_50s = (P.p30+P.p50+P.p80)/3.;
    P_10s = (P.p6+P.p8+P.p10+P.p13+P.p17)/5.;
    P_3s = (P.p2p2+P.p3+P.p4)/3.;
    P_1s = (P.p0p7+P.p1+P.p1p5)/3.;
    Pst = sqrt(0.0314*P.p0p1+0.0525*P_1s+0.0657*P_3s+0.28*P_10s+0.08*P_50s);
end