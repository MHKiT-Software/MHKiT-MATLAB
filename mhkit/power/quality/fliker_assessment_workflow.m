%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Conduct flicker assessment according to IECTS-30
%   Workflow of MATLAB functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. prepare data needed: 
% u_m, i_m, Sr, Un, In, SCR, fg, & fs
%
Sr = 3e6; Un=12e3; In=144;SCR=20;fg=60;fs=50e3;
% % OBS data: 230627_2019_11_11_18_10_PowRaw.tdms 
% fs = 50e3; fg=60Hz;
data = readmatrix('..\OneDrive - NREL\MHKiT-MATLAB\power\data\230627_2019_11_11_18_10_PowRaw.csv');
% % Resolve the time jump in the first 20000 OBSs:
xtime = data(20001:end,7);
xum=data(20001:end,1);xim=data(20001:end,4);
u_m = struct();u_m.time = xtime - xtime(1); u_m.data = xum;
i_m = struct();i_m.time = xtime - xtime(1); i_m.data = xim;
% calculate Un and In from u_m and i_m
Un = rms(u_m.data); In = rms(i_m.data);
%%


%% 1. Construct the fictitious grid:
[Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);
%%

%% 2. Calculate ideal voltage (u0) from measured voltage (u_m)
% 2.1 Calculate the fundamental frequency and alpha_0 of u_m:
% method STFT:
%method = {'Window',rectwin(M),'OverlapLength',L,...
%    'FFTLength',128,'FrequencyRange','twosided'};
%[alpha0,freq] = calc_alpha0_freq(u_m,'stft',method); 
% method ZCD (recommended):
[alpha0,freq] = calc_alpha0_freq(u_m,'zcd',{});
% 2.2 Calculate the electrical angle (alpha_m):
alpha_m = calc_electrical_angle(freq,alpha0);
% 2.3 Calculate ideal voltage (u0):
u0 = calc_ideal_voltage(Un,alpha_m);
%% check u0 to make sure it fullfills the requirements in the IECTS
hold off;
plot(u_m.time,u_m.data/max(u_m.data));hold on;
plot(u_m.time,u0/max(u0));legend('um','u0');xlim([0 1])

%% 3. Calculate simulated voltage (u_fic)
u_fic = calc_simulated_voltage(u0,i_m,Rfic,Lfic);
%%

%% 4. Calculate flicker emission (P_st,fic)
% P_stfic = calc_flicker_emission_fic(u_fic);
% prep input for digital flickermeter
% use digital flickermeter
P_stfic = Pst;
%%

%% 5. Calculate flicker coefficient
Xfic = 2*pi*fg*Lfic;
S_kfic = (Un^2)./sqrt(Rfic.^2+Xfic.^2);
c = calc_flicker_coefficient(P_stfic,S_kfic,Sr);