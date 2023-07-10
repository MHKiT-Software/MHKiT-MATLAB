%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduct flicker assessment according to IECTS-30
%   Workflow of MATLAB functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. prepare data needed: 
% u_m, i_m, Sr, Un, In, SCR, fg, & fs
Sr = 3e6; Un=12e3; In=144; SCR=20; fg=60; fs=50e3;
%% opt1: generate test data:
% IEC 61400-21-1: B.3.5 slow freq change
% DeltaI_I = [0 0 0 0];% pure sine wave
fm = 25;% modulating frequency 
% TableB.3,fg=60, SCR=50
% DeltaI_I = [8.040 9.376 11.479 12.844]; % fm=33.3
% TableB.2,fg=60, SCR=20
DeltaI_I = [4.763 5.726 7.640 9.488]; % fm=25
% DeltaI_I = [3.212 3.958 5.644 7.711]; % fm=20
[i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,DeltaI_I,5,650);
%% opt 2: OBS data: 230627_2019_11_11_18_10_PowRaw.tdms 
% % fs = 50e3; fg=60Hz;
% data = readmatrix('..\OneDrive - NREL\MHKiT-MATLAB\power\data\230627_2019_11_11_18_10_PowRaw.csv');
% % % Resolve the time jump in the first 20000 OBSs:
% xtime = data(20001:end,7);
% xum=data(20001:end,1);xim=data(20001:end,4);
% u_m = struct();u_m.time = xtime - xtime(1); u_m.data = xum;
% i_m = struct();i_m.time = xtime - xtime(1); i_m.data = xim;
% % calculate Un and In from u_m and i_m
% Un = rms(u_m.data); In = rms(i_m.data);
%%


%% 1. Construct the fictitious grid:
[Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
%%

%% 2. Calculate ideal voltage (u0) from measured voltage (u_m)
% 2.1 Calculate the fundamental frequency and alpha_0 of u_m:
% 2.1.1 method STFT:
% method = {'Window',rectwin(M),'OverlapLength',L,...
%    'FFTLength',128,'FrequencyRange','twosided'};
% opts = {'Window',rectwin(50000),'OverlapLength',25000,...
% 'FFTLength',50e3,'FrequencyRange','onesided'};
% [alpha0,freq] = calc_alpha0_freq(u_m,'stft',opts); 
% 2.1.2 method ZCD:
[alpha0,freq] = calc_alpha0_freq(u_m,'zcd',{});
% 2.2 Calculate the electrical angle (alpha_m):
alpha_m = calc_electrical_angle(freq,alpha0);
% 2.3 Calculate ideal voltage (u0):
u0 = calc_ideal_voltage(Un,alpha_m);
%% check u0 to make sure it fullfills the requirements in the IECTS
hold off;
plot(u_m.time,u_m.data);hold on; grid on;
plot(u_m.time,u0);legend('um','u0');xlim([10,10.5])

%% 3. Calculate simulated voltage (u_fic)
u_fic = calc_simulated_voltage(u0,i_m,Rfic,Lfic);
%%

%% 4. Calculate flicker emission (P_st,fic)
% run power_KALFlicker to get x_in, enable_in, & V_flicker_in
% prep input for flickermeter: u_fic_in
% alternative: from OBS data get u_fic, time from readmatrix()
size(u_fic)
time = u_m.time;
% check input to the digital flicker:
u_fic_in = timeseries(u_fic(:,2),time);
% run digital flickermeter in simulink
%% calculate Pst
S5 = out.S5;
% open statistical analyzer for flickermeter
%% 5. Calculate flicker coefficient
P_stfic = 0.0405;
Xfic = 2*pi*fg*Lfic;
S_kfic = (Un^2)./sqrt(Rfic.^2+Xfic.^2);
coef_flicker = calc_flicker_coefficient(P_stfic,S_kfic,Sr);
%% 6. Keep test data:
writematrix(i_m.data(1:100,:),"testdata/B.3.5_fm25_SCR20_im.txt");
writematrix(u_m.data(1:100),  "testdata/B.3.5_fm25_SCR20_um.txt");
%writematrix(Rfic,             "testdata/FicGrid_SCR50_fg60_Rfic.txt");
%writematrix(Lfic,             "testdata/FicGrid_SCR50_fg60_Lfic.txt");
writematrix(alpha0,           "testdata/B.3.5_fm25_SCR20_alpha0.txt");
writematrix(freq.data(1:100), "testdata/B.3.5_fm25_SCR20_freq.txt");
writematrix(alpha_m(1:100),   "testdata/B.3.5_fm25_SCR20_alpham.txt");
writematrix(u0(1:100),        "testdata/B.3.5_fm25_SCR20_u0.txt");
writematrix(u_fic(1:100,:),   "testdata/B.3.5_fm25_SCR20_ufic.txt");

%%

