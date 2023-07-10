function [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,opt,T)
%%%%%%%%%%%%%%%
% Generate test data according to IEC61400-21-1 Annex B.
% B.3.5 Slow frequency changes opt=5
% B.3.1 - B.3.4 other tests opt=2-4
% opt = 0: pure sine wave (alpha0=pi/6) for both i_m and u_m
%%%%%%%%%%%%%%%
    % time: up to 650s
    t = 0:1/fs:T-1/fs; %600=10min duration
    t = t';
    %% i_m (see TableB.2 for DeltaI/I)
    % SCR=20, fm = 20 
    i_m = struct();i_m.time = t;
    if opt ==0
        i_m.data = sqrt(2)*In*(1+1/100/2*sin(2*pi*fm*t)*DeltaI_I).*...
            sin(2*pi*fg*t);
        u_m = struct();u_m.time = t;
        u_m.data = sqrt(2/3)*Un*sin(2*pi*fg*t+pi/6);
        return
    end
    %% i_m
    i_opts = [2 3 4];
    disp(opt);
    if any(i_opts == opt)
        % Eq (B.1)
        i_m.data = sqrt(2)*In*(1+1/100/2*sin(2*pi*fm*t)*DeltaI_I).*...
            sin(2*pi*fg*t);
    elseif opt == 5
        % B.3.5 slow frequency changes: 
        % Eq (B.5) & Eq (B.6)
        func = fg+0.05.*sin(2*pi/60.*t);
        % plot(t,func);grid;
        i_m.data = sqrt(2)*In*(1+1/100/2*sin(2*pi*fm*t)*DeltaI_I).*...
            sin(2*pi*cumtrapz(t,func));
        %plot(t,func);xlim([0 T]);
    else
        ME = MException('MATLAB:gen_test_data',['option ' ...
            'does not exist']);
        throw(ME);
    end
    %% u_m
    u_m = struct();u_m.time = t;
    if opt == 2
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
        % NOTE: Eq(B.4) is '...(1+3+1/100*1/2*sin...'
        % fv: modulating frequency for voltage
        % The test is divided into 60 different cases, with fν in [0.5:0.5:30].
        fv = 0.5; 
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
