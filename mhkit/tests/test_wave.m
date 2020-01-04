%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This progam is written to test the wave module of MHKit Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import testing data from Wave Energy Prize

file='../../../../mhkit-data/Wave_Energy_Prize/iws1.txt';
T=readtable(file);
omega = 0.1:0.01:3.5;
test.f=omega/(2*pi);
test.Hs = 2.5;
test.Tp = 8;

testCase = matlab.unittest.TestCase.forInteractiveUse;

%testing elevation spectrum 
sample_rate = 50;
nnft = 1024;
Tarray=T(:,2:end);
Time_array=T(:,1);
wave_spectra=elevation_spectrum(Tarray,sample_rate,nnft,Time_array);
s1=size(wave_spectra.spectrum);
s2=size(Tarray);
assertEqual(testCase,s1(2),s2(2));

% Testing optional input parameters
%wave_spectra=elevation_spectrum(Tarray,sample_rate,nnft,Time_array,false);
%wave_spectra=elevation_spectrum(Tarray,sample_rate,nnft,Time_array,false,"parzen");
%begin intentional fails to test asserts
%wave_spectra=elevation_spectrum(Tarray,sample_rate,nnft,Time_array,false,"parzen",1);
%wave_spectra=elevation_spectrum(Tarray,sample_rate,nnft,Time_array,false,5);

%testing create_spectra.m 

pm=create_spectra('pierson_moskowitz_spectrum',test.f,test.Tp);
Tp0=peak_period(pm);
assertEqual(testCase,Tp0,test.Tp,'AbsTol',0.1);

pm=create_spectra('bretschneider_spectrum',test.f,test.Tp,test.Hs);
Hm0=significant_wave_height(pm);
assertEqual(testCase,Hm0.values,test.Hs,'AbsTol',0.1);

pm=create_spectra('jonswap_spectrum',test.f,test.Tp,test.Hs);
pm=create_spectra('jonswap_spectrum',test.f,test.Tp,test.Hs,2.4);

results= [0.1261;0.1386;0.2804;0.1073;0.3068;0.1551];

Hm0=significant_wave_height(wave_spectra);





