% 1. Construct the fititious grid:
[Rfic,Lfic]=calc_Rfic_Lfic(Sr,SCR,Un,fg);

% 2. Calculate ideal voltage (u0) from measured voltage (u_m)
% Calculate the fundamental frequency of u_m:
method = {'Window',rectwin(M),'OverlapLength',L,...
    'FFTLength',128,'FrequencyRange','twosided'};
[alpha0,freq] = calc_frequency(u_m,method); 
% Interpolate the frequencies to the time steps of u_m 
ft = struct();ft.time = u_m.time;
ft.data =interp1(seconds(freq.time),freq.data,...
    u_m.time,'linear','extrap');
% Calculate the ideal voltage:
u0 = calc_ideal_voltage(Un,alpha0,ft);

% 3. Calculate simulated voltage (u_fic)
u_fic = calc_simulated_voltage(u0,i_m,Rfic,Lfic);

% 4. Calculate flicker emission (P_st,fic)
P_stfic = calc_flicker_emission_fic(u_fic);

% 5. Calculate flicker coefficient
S_kfic = Un^2./sqrt(Rfic^2+Lfic^2);
c = calc_flicker_coefficient(P_stfic,S_kfic,Sr);