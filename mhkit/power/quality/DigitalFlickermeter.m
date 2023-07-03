% prep input for flickermeter
% from OBS data 
% u_fic = readmatrix("230630_ufic.txt");
% time = readmatrix("230630_time.txt");
% use gen_test_data to test B.3.5 from IECTS61400-21-1
size(u_fic)
time = u_m.time;
% check input to the digital flicker:
% run power_KALFlicker to get x_in, enable_in, & V_flicker_in
%%
u_fic_in = timeseries(u_fic(:,1),time);
% try 0.0001 s Ts
% u_fic_in = timeseries(u_fic(1:5:end,1),time(1:5:end));
% plot(u_fic_in)

%% digital flickermeter output P_st,fic
Pstfic = 0.1553;
Skfic = (Un^2)./sqrt(Rfic.^2+(2*pi*fg*Lfic).^2);
disp(Skfic/Sr);
coef_flicker = 20*Pstfic;
disp((coef_flicker-2)/2);