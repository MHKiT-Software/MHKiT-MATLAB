%% Example: MHKiT-MATLAB Power Module
%
% This example will demonstrate the <https://mhkit-software.github.io/MHKiT/mhkit-matlab/api.power.html
% MHKiT Power Module> functionality to compute power, instantaneous frequency,
% and harmonics from time series of voltage and current.

%% Load Power Data
% We will begin by reading in time-series data of measured three-phase (A, B,
% and C) voltage and current. IEC TS 62600-30 requires power quality assessments
% to use time series of at least 10 minutes, but for this example we will only
% look at a fraction of a second of data.

% Read in time-series data of voltage (V) and current (I)
power_table = readtable('./data/power/2020224_181521_PowRaw.csv');

% Display without scientific notation, then restore the default format
format shortG
power_table
format short

%%
% To use the MHKiT-MATLAB power module we need to create structures of current
% and voltage.

current.time = posixtime(power_table.Time_UTC); % setting time in current structure
voltage.time = posixtime(power_table.Time_UTC); % setting time in voltage structure
T2 = mergevars(power_table, [2 3 4]); % combining the voltage time series into one table variable
T3 = mergevars(T2, [3 4 5]); % combining the current time series into one table variable
current.current = T3.Var3;
voltage.voltage = T3.Var2;

%%
% IEC TS 62600-30 clause 6.1 requires the generator sign convention: power
% flowing from the marine energy converter to the grid is positive. In this
% dataset every sample of the instantaneous power (the sum of voltage times
% current over the three phases) is negative, which indicates the current
% transformers were connected with the opposite orientation. We therefore
% flip the sign of the measured currents so that the exported power is positive.

current.current = -current.current;

%%
% Before computing anything, let's look at the first 0.1 s of the three-phase
% voltage and current signals.

% Time since the start of the record [s]
elapsed_time = voltage.time - voltage.time(1);
first_tenth = elapsed_time <= 0.1;

figure('Position', [100, 100, 1600, 600]);
plot(elapsed_time(first_tenth), voltage.voltage(first_tenth, :));
title('Three-Phase Voltage');
xlabel('Time [s]');
ylabel('Voltage [V]');
ax = gca;
ax.YAxis.Exponent = 0; % show plain volts instead of scientific notation
legend('Phase A', 'Phase B', 'Phase C');
grid on;

figure('Position', [100, 100, 1600, 600]);
plot(elapsed_time(first_tenth), current.current(first_tenth, :));
title('Three-Phase Current');
xlabel('Time [s]');
ylabel('Current [A]');
legend('Phase A', 'Phase B', 'Phase C');
grid on;

%% Power Characteristics
%
% The MHKiT |power/characteristics| submodule can be used to compute basic quantities
% of interest from voltage and current time series. In this example we will calculate
% active AC power and fundamental frequency, using the loaded voltage and current
% time-series.
%
% To compute the active AC power, we will need the power factor, the ratio of
% active (real) power to apparent power ($PF = P/S$, between 0 and 1). It
% describes how effectively the voltage and current at the connection deliver
% real power; it is not the device's conversion efficiency ($P_{out}/P_{in}$).
% For a grid-connected marine energy converter, the power factor at the grid
% connection is largely set by the power electronics' reactive power control
% (for example constant power factor, constant reactive power, or voltage
% control modes; IEC TS 62600-30 clause 6.7 and Table A.5). Here we assume a
% constant value of 0.96.

% Set the power factor for the system
power_factor = 0.96;

% Compute the instantaneous AC power in watts
ac_power = ac_power_three_phase(voltage, current, power_factor);

% Display the result in kilowatts
figure('Position', [100, 100, 1600, 600]);
plot(datetime(ac_power.time, "ConvertFrom", "posixtime"), ac_power.power / 1000);
title('AC Power');
ylabel('Power [kW]');
xlabel('Time');
ax = gca;
ax.YAxis.Exponent = 0; % show plain kilowatts instead of scientific notation

%%
% *Fundamental Frequency*
%
% Using the 3 phase voltage measurements we can compute the fundamental frequency
% of the voltage time series. The time-varying fundamental frequency is a required
% metric for power quality assessments. The function |calc_fundamental_freq()|
% provides two methods for calculating the fundamental frequency:
%
% # Short-Time Fourier Transform (STFT)
% # Zero-Crossing Detection (ZCD)
%
% Here we illustrate both.
%
% First, let's look at the STFT method. To obtain an accurate result the STFT
% window settings are crucial, so some exploratory analysis with different sets
% of parameters may be needed.
%
% Note: The STFT method requires the Signal Processing Toolbox.

% Hilbert method: instantaneous_frequency
% inst_freq = instantaneous_frequency(voltage)

% STFT method:
% Check that the Signal Processing Toolbox functions are installed. license('test')
% only checks that a license exists, not that the toolbox is installed.
has_signal_toolbox = exist('stft', 'file') ~= 0 && exist('rectwin', 'file') ~= 0;
if has_signal_toolbox
    % Set up method options
    methodopts = {
        'Window', int32(5000), ...
        'OverlapLength', 3750, ...
        'FFTLength', 50e3, ...
        'FrequencyRange', 'onesided'
    };

    % Prep input and output
    u_m = struct();
    u_m.time = voltage.time;
    fund_freq = struct();
    fund_freq.time = voltage.time;
    fund_freq.data = zeros(size(voltage.voltage));

    for i = 1:3
        u_m.data = voltage.voltage(:,i);
        [~, freq] = calc_fundamental_freq(u_m, 'stft', methodopts);
        fund_freq.data(:,i) = freq.data;
    end
    fund_freq

    % Display the result
    figure('Position', [100, 100, 1600, 600]);
    plot(datetime(fund_freq.time, "ConvertFrom", "posixtime"), fund_freq.data);
    title('Fundamental Frequency: STFT');
    ylabel('Frequency [Hz]');
    xlabel('Time');
    ylim([50, 70]);
else
    warning('Signal Processing Toolbox is not installed. Skipping STFT analysis...');
end

%%
% Now, let's look at the ZCD method. Note that the ZCD method may not work for
% distorted voltage with multiple zero crossings, such as the case described in
% IEC 61400-21-1:2019 Annex B.3.3. As stated above, exploratory analysis is needed.

% ZCD method:
% Prep input and output
u_m = struct();
u_m.time = voltage.time;
fund_freq = struct();
fund_freq.time = voltage.time;
fund_freq.data = zeros(size(voltage.voltage));

% For the zcd method, we don't need to worry about the method options
methodopts = {};

for i = 1:3
    u_m.data = voltage.voltage(:,i);
    [~, freq] = calc_fundamental_freq(u_m, 'zcd', methodopts);
    fund_freq.data(:,i) = freq.data;
end
fund_freq

figure('Position', [100, 100, 1600, 600]);
plot(datetime(fund_freq.time, "ConvertFrom", "posixtime"), fund_freq.data);
title('Fundamental Frequency: ZCD');
ylabel('Frequency [Hz]');
xlabel('Time');
ylim([50, 70]);

%% Power Quality
%
% The power quality submodule can be used to compute voltage fluctuations
% (flicker), and harmonics and harmonic distortion of current and voltage,
% following IEC TS 62600-30 and IEC 61000-4-7. Harmonics and harmonic distortion
% are required as part of a power quality assessment and characterize the
% quality of the power being produced. We start with the current harmonics.

% Set the sampling frequency of the dataset
sample_freq = 50000; % [Hz]

% Set the frequency of the grid the device would be connected to
grid_freq = 60; % [Hz]

% Set the rated current of the device
rated_current = 18.8; % [A]

% Calculate the current harmonics
h = harmonics(current, sample_freq, grid_freq);

% Display the results
figure('Position', [100, 100, 1600, 600]);
plot(h.harmonic, h.amplitude);
title('Current Harmonics');
xlabel('Frequency [Hz]');
ylabel('Harmonic Amplitude [A]');
legend('Phase A', 'Phase B', 'Phase C');
xlim([0, 900]);

%%
% *Harmonic Subgroups*
%
% IEC TS 62600-30 clause 7.3 requires the harmonic currents up to 50 times the
% grid frequency to be reported as harmonic subgroups. We calculate them from
% the harmonics and the grid frequency.

% Calculate harmonic subgroups
h_s = harmonic_subgroups(h, grid_freq);

% Display the first ten harmonic orders (row 1 is DC, row 2 is the fundamental)
harmonic_order = (0:9)';
array2table([harmonic_order, h_s.harmonic(1:10), h_s.amplitude(1:10, :)], ...
    'VariableNames', {'Order', 'Frequency [Hz]', 'Phase A [A]', 'Phase B [A]', 'Phase C [A]'})

%%
% *Total Harmonic Current Distortion (THCD)*
%
% Finally we compute the total harmonic current distortion as a percentage.
% Matching MHKiT-Python, the subgrouped harmonics of order 2 to 49 are summed
% and normalized by the fundamental subgroup, and the result is for the first
% column (Phase A) only. The rated current is accepted but not yet used; the
% IEC TS 62600-30 formula (7) normalization by rated current and per-phase
% output are under review.

THCD = total_harmonic_current_distortion(h_s, rated_current)

%%
% *Flicker Assessment*
%
% Calculate the flicker coefficient following the steps in IEC TS 62600-30 section
% 7.2.1: MV connected marine energy converters.
%
% The first step is to calculate the simulated voltage of the fictitious grid
% (u_fic). There are two options for calculating u_fic.
%
% First, the user can use |flicker_ufic_workflow()|, which wraps up all the steps
% needed to calculate u_fic from the measured voltage (u_m) and current (i_m).
% All intermediate values are also returned to help with debugging.
%
% Alternatively, the user can call the corresponding functions sequentially. For
% illustration, we show both ways.
%
% Firstly, let's use |flicker_ufic_workflow()|.

% Prep the input: look at the first phase for an example
u_m = struct();
u_m.time = voltage.time;
u_m.data = voltage.voltage(:,1);
i_m = struct();
i_m.time = current.time;
i_m.data = current.current(:,1);

Sr = 4.2e5; % rated apparent power [VA]
Un = 14e3; % RMS value of the nominal voltage [V]
SCR = 20; % short-circuit ratio
fg = 60; % nominal grid frequency [Hz]

% Method settings for calculation of fundamental frequency
% of u_m, see above
method = 'zcd';
methodopts = {};

out = flicker_ufic_workflow(Sr, Un, SCR, fg, ...
    u_m, i_m, method, methodopts);
u0 = out.u0;
out

%%
% Calculating u_fic can be tricky because of the ideal voltage source (u0). As
% stated in IEC TS 62600-30 clause 7.2.2, u0 should fulfill the following two
% requirements:
%
% (1) be without any fluctuations
% (2) have the same electrical angle as the fundamental of u_m
%
% Therefore, we compare u0 and u_m to check that our output fulfills these
% requirements.

% Check u0
figure('Position', [100, 100, 1600, 600]);
plot(datetime(u_m.time, "ConvertFrom", "posixtime"), u_m.data);
hold on;
plot(datetime(u_m.time, "ConvertFrom", "posixtime"), u0);
hold off;
legend('u_m', 'u_0');
xlabel('Time');
ylabel('Voltage [V]');

%%
% Secondly, we can also calculate u_fic step by step. The calculation of u_fic
% can be achieved sequentially in 5 steps:

%%
% Step 1. Construct the fictitious grid: calculate resistance and inductance
% using |calc_Rfic_Lfic()|.

[Rfic, Lfic] = calc_Rfic_Lfic(Sr, SCR, Un, fg);

%%
% Step 2. Calculate the fundamental frequency and alpha_0:
%   -opt1: method ZCD: method = 'ZCD'; methodopts={}
%   -opt2: method STFT: method = 'stft'; methodopts = {'Window',rectwin(M),...
%           'OverlapLength',L,'FFTLength',128,'FrequencyRange','onesided'}

if has_signal_toolbox
    method = 'stft';
    methodopts = {
        'Window', rectwin(5000), ...
        'OverlapLength', 3750, ...
        'FFTLength', 50e3, ...
        'FrequencyRange', 'onesided'
    };
else
    warning('Signal Processing Toolbox is not installed. Using ZCD method...');
    method = 'ZCD';
    methodopts = {};
end

[alpha0, freq] = calc_fundamental_freq(u_m, method, methodopts);

%%
% Step 3. Calculate the electrical angle (alpha_m) of the fundamental of u_m
% using |calc_electrical_angle()|.

alpha_m = calc_electrical_angle(freq, alpha0);

%%
% Step 4. Calculate u0 from alpha_m and nominal voltage (Un) using |calc_ideal_voltage()|.

u0 = calc_ideal_voltage(Un, alpha_m);

%%
% Step 5. Calculate u_fic using |calc_simulated_voltage()|.

u_fic = calc_simulated_voltage(u0, i_m, Rfic, Lfic);

%%
% Again, we check our derived u0 to see if it fulfills the requirements.

% Check u0
figure('Position', [100, 100, 1600, 600]);
plot(datetime(u_m.time, "ConvertFrom", "posixtime"), u_m.data);
hold on;
plot(datetime(u_m.time, "ConvertFrom", "posixtime"), u0);
hold off;
legend('u_m', 'u_0');
xlabel('Time');
ylabel('Voltage [V]');

%%
% *Calculate the flicker emission value (P_stfic) from u_fic*
%
% Input u_fic into an appropriate digital flickermeter to get P_stfic.
%
% According to the standard, the evaluation time is 10 min, so the data must be
% at least 10 min long. Depending on the performance of the chosen flickermeter,
% the first 20 s may need to be discarded.
%
% u_fic has one column for each fictitious grid impedance phase angle, psi_k =
% 30, 50, 70, and 85 degrees (see |calc_Rfic_Lfic()|). IEC TS 62600-30 requires
% the flicker coefficient to be reported for all four angles (Annex A, Table
% A.6); here we illustrate a single angle, psi_k = 50 degrees.

% Select the fictitious grid impedance phase angle
psi_k = [30, 50, 70, 85]; % [degrees], column order of u_fic, Rfic, and Lfic
psi_idx = find(psi_k == 50);

% Prep input for the digital flickermeter
time = u_m.time - u_m.time(1);
u_fic_in = timeseries(u_fic(:, psi_idx), time);
u_fic_in

%%
% A digital flickermeter that also contains an interface of statistical analysis
% in simulink can be used to generate the instantaneous flicker level (|out.S5|)
% and then calculate the Pst using the statistical tool associated with it. More
% details can be found <https://www.mathworks.com/help/sps/powersys/ref/digitalflickermeter.html
% here>.
%
% Note that the user needs to specify RMS value for voltage, sampling frequencies,
% and other parameters to achieve a valid result. In addition, there may be some
% unreliable peaks in the instantaneous flicker levels derived from this particular
% flickermeter, when performing statistical analysis, be sure to discard the S5
% data at the very beginning of it.

% Prep input for the statistical tool
% out.S5 is the result generated from the digital flickermeter.
% S5 = out.S5;

%%
% After getting the instantaneous flicker levels, double-click the digital flickermeter
% to open statistical tools to calculate P_stfic.

%%
% *Calculate the flicker coefficient*
%
% Provide the flickermeter output (P_stfic) for the selected impedance phase
% angle. In this example the data duration is less than 1 s, so we use a
% placeholder value for P_stfic.

P_stfic = 0.10016;
Xfic = 2*pi*fg*Lfic;
S_kfic = (Un^2)./sqrt(Rfic.^2 + Xfic.^2); % IEC TS 62600-30 formula (5)
coef_flicker = calc_flicker_coefficient(P_stfic, S_kfic(psi_idx), Sr)

%%
% *Test the performance of the flicker assessment process*
%
% This part of the example illustrates how to test the performance of the user's
% flicker assessment process, including methods to generate simulated voltage
% of the fictitious grid and the performance of the (digital) flickermeter, following
% the guidance provided in IEC 61400-21-1:2019 Annex B.3: Verification test of
% the measurement procedure for flicker.
%
% The verification is straightforward: we generate test data with predetermined
% flicker coefficients and verify the flicker assessment process by comparing
% the derived flicker coefficients with the predetermined ones. The generated
% test data are measured voltage and current time series (|u_m| and |i_m|). The
% user can use |gen_test_data()| to generate test data for five scenarios: (1) a
% pure sine wave (coefficient = 0) and (2)-(5) the scenarios in Annex B.3.2 to
% B.3.5, and test them one by one. Here we illustrate |gen_test_data()| by
% generating a distorted u_m(t) with multiple zero crossings (Annex B.3.3).

% Set up proper input parameters:
% scenario 3 with duration of 10s
opt = 3; T=10;
% rated parameters: voltage [V], current [A], sample rate [Hz], apparent power [VA]
Un=12e3; In=144; fs=50e3; Sr = 3e6;
% other parameters for generating distorted data
fg=60;  SCR=20; fm=25;
fv=0.5; % fv only needed for scenario B.3.4
% According to IEC 61400-21-1:2019 Table B.2:
DeltaI_I = [4.763 5.726 7.640 9.488];
[i_m,u_m]=gen_test_data(Un,In,fg,fs,fm,fv,DeltaI_I,opt,T)

%%
% The generated i_m has four time series (columns), one for each impedance
% phase angle psi_k = 30, 50, 70, and 85 degrees. u_m is the same for all four
% angles, so it has a single time series.

figure('Position', [100, 100, 1600, 600]);
plot(i_m.time, i_m.data(:,1));
xlim([0 1]);
xlabel("Time [s]");
ylabel("i_m [A]");
title("Generated current, \psi_k = 30 degrees");

figure('Position', [100, 100, 1600, 600]);
plot(u_m.time, u_m.data);
xlim([0 0.1]);
xlabel("Time [s]");
ylabel("u_m [V]");
title("Generated voltage with multiple zero crossings");
grid on;
