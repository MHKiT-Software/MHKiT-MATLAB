%%% function
clear; close all;

%% import params

fpath = '../../../examples/data/loads/mler_data.csv';

wave_freq = linspace(0,1,500);
wave_spectrum = jonswap_spectrum(wave_freq,15.1,9);
validation = readtable(fpath);

RAO = validation.RAO;

%% processing

% convert from Hz to rad/s
freq = wave_spectrum.frequency * (2*pi);
wave_spec = wave_spectrum.spectrum / (2*pi);
dw = (2*pi - 0) / (length(freq)-1);

% allocate variables
S_R = zeros(length(freq));
zS = zeros(length(freq));
zA = zeros(length(freq));
zCoeffA_Rn = zeros(length(freq));
zphase = zeros(length(freq));

% Note: waves.A is "S" in Quon2016; 'waves' naming convention matches WEC-Sim conventions (EWQ)
% Response spectrum [(response units)^2-s/rad] -- Quon2016 Eqn. 3


%% 

S_R(:) = abs(RAO)^2 * (2*wave_spec);



