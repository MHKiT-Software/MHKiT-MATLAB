function [i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I)
%%%%%%%%%%%%%%%
% Generate test data according to IEC61400-21-1 Annex B.
% B.3.5 Slow frequency changes
%
%%%%%%%%%%%%%%%
    % B.3.5 Slow frequency changes
    % Sr = 3e6; Un = 12e3; In = 144;
    % fs = 50e3; fg = 60; 
    t = 0:1/fs:600-1/fs; %600=10min duration
    t = t';
    %fm = 20;%modulating frequency (TableB.2,fg=60)
   
    func = fg+0.05.*sin(2*pi/60.*t);
    %%
    plot(t,func);grid;
    %% u_m
    u_m = struct();u_m.time = t;
    u_m.data = sqrt(2/3)*Un*(1+0.25/100/2*sin(2*pi*8.8*t)).*...
        sin(2*pi*cumtrapz(t,func));
    %%
    % plot(t,u_m.data);xlim([0 1]);grid;
    %% i_m
    % TableB.2,fg=60, SCR=20, fm=20
    % SCR=20;fm = 20;%modulating frequency 
    % DeltaI_I = 3.212;% for Phik=30
    i_m = struct();i_m.time = t;
    i_m.data = sqrt(2)*In*(1+DeltaI_I/100/2*sin(2*pi*fm*t)).*...
        sin(2*pi*cumtrapz(t,func));


    % briefly check:
    %   [i_m1,u_m1] = gen_test_data(12e3,144,60,50e3,20,3.212);
    %   isequal(i_m1,i_m)
    %   
    %   ans =
    %   
    %     logical
    %   
    %      1
    %   
    %   isequal(u_m1,u_m)
    %   
    %   ans =
    %   
    %     logical
    %   
    %      1

end
