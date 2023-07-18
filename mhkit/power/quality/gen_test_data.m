function [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Generate test data according to IEC61400-21-1 Annex B.3 to be used for
%   the verification test of the measurement procedure for flicker.
% 
% Parameters
% -----------
%   Un: double 
%       RMS value of the nomianal voltage of the grid (V).
%   In: double 
%       Nominal current (A).
%   fg: double 
%       Nominal grid frequency, 60Hz or 50Hz
%   fs: double 
%       Sampling frequency (Hz).
%   fm: double 
%       Modulating frequency for current.
%   fv: double 
%       Modulating frequency for voltage used in B.3.4. The test is divided 
%       into 60 different cases, with fν in [0.5:0.5:30].
%   DeltaI_I: double array (4) 
%       Relative current changes (%) according to Table B.2 and Table B.3
%       in IECTS61400-21-1.
%   opt: int 
%       Option number used to indicate different test datasets. 
%   T: double 
%       Duration (s) of the generated u_m and i_m.
%  
% Returns
% -------
%   u_m: struct() 
%       u_m(t) is the measured instantaneous voltage (V)
%       .time: time at each measurement (s)
%       .data: array of size (ntime) measured instantaneous voltage (V).
%   i_m: struct()
%       i_m(t) is the measured instantaneous current (A)
%       .time: time at each measurement
%       .data: array of size (ntime, 4) measured instantaneous current (A)
%       at impedance phase angles = 30, 50, 70, 85.
% Note
% -------
% 1. Options (opt=0,[2,5])
% opt = 0: Pure sine wave i_m and u_m with alpha0=pi/6.
% opt = 2 - 5
% B.3.2 Fictitious grid performance testing, pure sine wave for u_m
% B.3.3 Distorted u_m with multiple zero crossings
% B.3.4 Distorted u_m with inter-harmonic modulation 
% B.3.5 Slow frequency changes
%%%%%%%%%%%%%%%
    %% time: 
    t = 0:1/fs:T-1/fs; %600=10min
    t = t';
    %% i_m (see TableB.2 & TableB.3 DeltaI/I)
    i_m = struct();i_m.time = t;
    i_opts = [0 2 3 4];
    % disp(opt);
    if any(i_opts == opt)
        % Eq (B.1)
        i_m.data = sqrt(2)*In*(1+1/100/2*sin(2*pi*fm*t)*DeltaI_I).*...
            sin(2*pi*fg*t);
    elseif opt == 5
        % B.3.5 slow frequency changes: 
        % Eq (B.5) & Eq (B.6)
        func = fg+0.05.*sin(2*pi/60.*t);
        i_m.data = sqrt(2)*In*(1+1/100/2*sin(2*pi*fm*t)*DeltaI_I).*...
            sin(2*pi*cumtrapz(t,func));
        %plot(t,func);xlim([0 T]);grid on;
    else
        ME = MException('MATLAB:gen_test_data',['option ' ...
            'does not exist']);
        throw(ME);
    end
    %% u_m
    u_m = struct();u_m.time = t;
    if opt ==0
        % Pure sine wave
        u_m.data = sqrt(2/3)*Un*sin(2*pi*fg*t+pi/6);
    elseif opt == 2
        % B.3.2 fictitious grid performace testing: Eq (B.2)
        u_m.data = sqrt(2/3)*Un*sin(2*pi*fg*t);
    elseif opt ==3
        % B.3.3 distorted u_m(t) voltage with multiple zero crossings: Eq (B.3)
        % Table B.4:
        percents = [5 6 5 1.5 3.5 3 2 1.76 1.41 1.27 1.06 0.97];
        horder = [3 5 7 9 11 13 17 19 23 25 29 31];
        u_m.data = sqrt(2/3)*Un*(1+0.25/100/2*sin(2*pi*8.8*t)).*...
            (sin(2*pi*fg*t)+sum(percents/100.*sin(2*pi*fg*t.*horder+pi),2));
    elseif opt==4    
        % B.3.4 distorted u_m(t) voltage with inter-harmonic modulation: Eq (B.4)
        % NOTE: Eq(B.4) is '...(1+3+1/100*1/2*sin...' in IECTS, maybe wrong. 
        u_m.data = sqrt(2/3)*Un*(1+3/100/2.*sin(2*pi*fv*t)).*sin(2*pi*fg*t);
    elseif opt==5
        % B.3.5 slow frequency changes: Eq (B.7)
        u_m.data = sqrt(2/3)*Un*(1+0.25/100/2*sin(2*pi*8.8*t)).*...
            sin(2*pi*cumtrapz(t,func));
        % plot(t,u_m.data);xlim([0 1]);grid;
    else
        ME = MException('MATLAB:gen_test_data',['option ' ...
            'does not exist']);
        throw(ME);
    end

end
