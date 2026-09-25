function v = spectral_width(S, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave spectral width from spectra
%
% Computed as sqrt(m0 * m_{-2} / m_{-1}^2 - 1) where m0, m_{-1}, and m_{-2}
% are spectral moments (Eq 16 in IEC 62600-101 Ed. 2.0 en 2024).
%
% Parameters
% ------------
% S : struct, numeric, table, or timetable
%   struct: S.spectrum (vector or matrix [m^2/Hz]), S.frequency (vector
%     [Hz]), optional S.time
%   numeric: spectral density array, frequency required as second argument
%   table: one variable per spectrum (rows), frequency required as second
%     argument, optional 'time' variable
%   timetable: one variable per spectrum (rows = RowTimes), frequency
%     required as second argument
% frequency_bins : vector (optional)
%   Bin widths for frequency of S. Required for unevenly sized bins.
%
% Returns
% ---------
% v : double, column vector, table, or timetable [dimensionless]
%   Spectral width, one value per spectrum. Matches the container style
%   of S: numeric for struct/numeric input, table/timetable for
%   table/timetable input.
%
% Examples
% --------
%     CDIP Example: Station Number 225: Kaneohe Bay, WETS, Oahu, HI
%     https://cdip.ucsd.edu/themes/cdip?pb=1&u2=s:225:st:1&d2=p70
%     >> station_number = '225';
%     >> data_type = 'realtime';
%     >> years = 2025;
%     >> parameters = {'waveEnergyDensity'};
%     >> data = cdip_request_parse_workflow('station_number', station_number, 'data_type', data_type, 'years', years, 'parameters', parameters);
%     >> frequency = data.metadata.wave.waveFrequency;  % [Hz], 64x1 single
%     frequency =
%         0.0250
%         0.0300
%         0.0350
%          :
%         0.5600
%         0.5700
%         0.5800
%
%     >> % waveEnergyDensity is stored [time x frequency]; transpose to
%     >> % match MHKiT's frequency-as-rows, spectra-as-columns convention.
%     >> spectrum = data.data.wave2D.waveEnergyDensity';  % [m^2/Hz], 64x17520 double
%     spectrum =
%        0.0002     0.0001     0.0002  ...     0.0006     0.0003     0.0004
%        0.0007     0.0003     0.0005  ...     0.0021     0.0007     0.0010
%        0.0032     0.0030     0.0020  ...     0.0047     0.0029     0.0038
%          :
%        0.0165     0.0122     0.0138  ...     0.0186     0.0079     0.0128
%        0.0124     0.0116     0.0080  ...     0.0128     0.0120     0.0088
%        0.0184     0.0097     0.0116  ...     0.0095     0.0084     0.0080
%
%     >> time = data.data.wave.waveTime;  % 17520x1 datetime
%     time =
%     01-Jan-2025 00:00:00
%     01-Jan-2025 00:30:00
%     01-Jan-2025 01:00:00
%          :
%     31-Dec-2025 22:30:00
%     31-Dec-2025 23:00:00
%     31-Dec-2025 23:30:00
%
%     Struct (CDIP real-world data): single (most recent) spectrum
%     >> S.frequency = frequency;
%     >> S.spectrum = spectrum(:,end);
%     >> v = spectral_width(S);
%     v =
%         0.4207  % [dimensionless]
%
%     Numeric (CDIP real-world data): matrix, one spectrum per column
%     >> v = spectral_width(spectrum, frequency);
%     v =  % [dimensionless]
%         0.3669
%         0.3937
%         0.3947
%          :
%         0.4155
%         0.3892
%         0.4207
%
%     Table (CDIP real-world data): one row per spectrum, one variable per frequency bin
%     >> col_names = mhkit_frequency_to_column_names(frequency);
%     >> T = array2table(spectrum', 'VariableNames', col_names);
%     >> v_table = spectral_width(T, frequency);
%     v_table =  % [dimensionless]
%     17520×1 table
%         spectral_width
%         ______________
%            0.36695    
%            0.39367    
%             0.3947    
%         :
%            0.41548    
%             0.3892    
%            0.42073
%
%     Timetable (CDIP real-world data): RowTimes carried through to the output
%     >> TT = array2timetable(spectrum', 'RowTimes', time, 'VariableNames', col_names);
%     >> v_tt = spectral_width(TT, frequency);
%     v_tt =  % [dimensionless]
%     17520×1 timetable
%                 time            spectral_width
%         ____________________    ______________
%         01-Jan-2025 00:00:00       0.36695    
%         01-Jan-2025 00:30:00       0.39367    
%         01-Jan-2025 01:00:00        0.3947    
%         :
%         31-Dec-2025 22:30:00       0.41548    
%         31-Dec-2025 23:00:00        0.3892    
%         31-Dec-2025 23:30:00       0.42073
%
%     WEC-Sim Output Example
%     >> S = load('examples/data/RM3MooringMatrix_matlabWorkspace.mat', 'output');
%     >> elevation = S.output.wave.elevation;  % [m], RM3 float, 40001x1 double
%     >> raw_time = S.output.wave.time;  % [s], 40001x1 double
%     >> sample_rate = 1 / (raw_time(2) - raw_time(1));  % [Hz], 100
%     >> % Note: IEC 62600-101 Ed. 2.0 en 2024, "Wave energy resource
%     >> % assessment and characterization", section 6.5.2 specifies a
%     >> % wave record length of minimum 1200 s (20 min), up to 3600 s
%     >> % (60 min) for better spectral resolution/precision. This
%     >> % WEC-Sim record is only ~6.7 min, so treat this as illustrative,
%     >> % not a fully representative sea-state estimate.
%     >> %
%     >> % Note: spectral statistics assume an irregular wave record.
%     >> % They are not meant for regular, single-frequency (sine) wave
%     >> % tests, which WEC-Sim is often used to run.
%     >> Sxx = elevation_spectrum(elevation, sample_rate, 1000, raw_time);
%     >> frequency = Sxx.frequency;  % [Hz], 501x1 double
%     frequency =
%         0.0000
%         0.1000
%         0.2000
%          :
%        49.8000
%        49.9000
%        50.0000
%
%     >> spectrum = Sxx.spectrum;  % [m^2/Hz], 501x1 double
%     spectrum =
%         0.1439
%         0.8485
%         0.8739
%          :
%         0.0000
%         0.0000
%         0.0000
%
%     Struct (WEC-Sim output): the single spectrum
%     >> S.frequency = frequency;
%     >> S.spectrum = spectrum(:,end);
%     >> v = spectral_width(S);
%     v =
%         0.4031  % [dimensionless]
%
%     Numeric (WEC-Sim output): matrix, one spectrum per column
%     >> v = spectral_width(spectrum, frequency);
%     v =
%         0.4031  % [dimensionless]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    S
end

arguments (Repeating)
    varargin
end

arguments (Output)
    v
end

[spectrum, frequency, time, input_style, remaining] = mhkit_standardize_spectrum_input(S, 'spectral_width', varargin{:});

% Eq 16 in IEC 62600-101 Ed. 2.0 en 2024
if isempty(remaining)
    m0 = frequency_moment(spectrum, 0, frequency);
    mn1 = frequency_moment(spectrum, -1, frequency);
    mn2 = frequency_moment(spectrum, -2, frequency);
else
    m0 = frequency_moment(spectrum, 0, frequency, remaining{1});
    mn1 = frequency_moment(spectrum, -1, frequency, remaining{1});
    mn2 = frequency_moment(spectrum, -2, frequency, remaining{1});
end

v = sqrt(((m0 .* mn2) ./ (mn1.^2)) - 1);
v = v(:);
mhkit_verify_is_column_vector(v, 'function_name', mfilename);

v = mhkit_restore_spectrum_output(v, input_style, 'spectral_width', time);

end
