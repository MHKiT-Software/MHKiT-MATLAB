function harmonics = harmonics(x, freq, grid_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the harmonics from time series of voltage or current based on IEC 61000-4-7.
%
% Parameters
% ------------
%   x: structure
%       x.current : Current time series data [A] 
%       x.voltage : Voltage time series data [V]
%       x.time : Time vector [s]
%
%   freq: double
%       Frequency of the time-series data [Hz]
%
%   grid_freq: double
%       Value indicating if the power supply is 50 or 60 Hz. Options = 50 or 60
%
% Returns
% ---------
%   harmonics: structure
%       harmonics.amplitude : Harmonic amplitude values
%       harmonics.harmonic : Harmonic frequency values [Hz]
%       harmonics.type : Type of signal analyzed ('current' or 'voltage')
%
% Key Equations
% -------------
% 1. Sample spacing calculation:
%    sample_spacing = 1 / freq
%
% 2. FFT amplitude calculation:
%    harmonics_amplitude = abs(fft(signal_data))
%
% 3. Normalization:
%    normalized_amplitude = harmonics_amplitude / length(signal_data) * 2
%
% 4. Frequency bin calculation:
%    frequency_bins = (0:length(signal_data)-1) * freq / length(signal_data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStruct = @(x) isstruct(x);
    validNumeric = @(x) isnumeric(x) && isscalar(x);
    validGridFreq = @(x) isnumeric(x) && isscalar(x) && (x == 50 || x == 60);
    
    % Add required parameters
    addRequired(p, 'x', validStruct);
    addRequired(p, 'freq', validNumeric);
    addRequired(p, 'grid_freq', validGridFreq);
    
    % Parse inputs
    parse(p, x, freq, grid_freq);
    
    % Extract validated inputs
    x = p.Results.x;
    freq = p.Results.freq;
    grid_freq = p.Results.grid_freq;
    
    % Validate input structure has required fields
    if ~isfield(x, 'time')
        error('MHKiT:harmonics: x structure must contain time field');
    end
    
    % Determine signal type and extract data
    if isfield(x, 'current')
        signal_data = x.current;
        signal_type = 'current';
    elseif isfield(x, 'voltage')
        signal_data = x.voltage;
        signal_type = 'voltage';
    else
        error('MHKiT:harmonics: x structure must contain either current or voltage field');
    end
    
    % Validate frequency parameter
    if freq <= 0
        error('MHKiT:harmonics: freq must be positive');
    end
    
    % Get data dimensions
    data_size = size(signal_data);
    data_length = data_size(1);
    num_columns = data_size(2);
    
    % Validate time vector dimensions
    if length(x.time) ~= data_length
        error('MHKiT:harmonics: time vector length must match signal data length');
    end
    
    % Calculate sample spacing
    sample_spacing = 1.0 / freq;
    
    % Calculate frequency bin centers (only positive frequencies)
    frequency_bin_centers = (0:data_length-1) * freq / data_length;
    
    % Calculate FFT amplitude for each column
    harmonics_amplitude = zeros(data_length, num_columns);
    
    for col = 1:num_columns
        % Calculate FFT
        fft_result = fft(signal_data(:, col));
        
        % Calculate amplitude
        harmonics_amplitude(:, col) = abs(fft_result);
    end
    
    % Calculate Nyquist frequency and corresponding bin
    nyquist_freq = freq / 2;
    nyquist_bin = floor(data_length / 2) + 1;
    
    % Define parameters based on IEC 61000-4-7
    max_harmonic_base = 51; % Base harmonic number from IEC 61000-4-7
    harmonic_step = 5; % Step size in Hz for harmonic bins

    if grid_freq == 60
        max_freq = 3060; % Exactly 51st harmonic (60 * 51)
    elseif grid_freq == 50
        max_freq = 2570; % 51.4th harmonic (intentional extension)
    end

    % Create frequency range (exclusive upper bound like Python's arange)
    hz_range = 0:harmonic_step:(max_freq - harmonic_step);
    
    % Interpolate to standard frequency grid using nearest neighbor
    harmonics_reindexed = zeros(length(hz_range), num_columns);
    
    for col = 1:num_columns
        for i = 1:length(hz_range)
            target_freq = hz_range(i);
            
            % Skip frequencies above Nyquist frequency to avoid aliasing
            if target_freq > nyquist_freq
                harmonics_reindexed(i, col) = 0;
                continue;
            end
            
            % Calculate the exact FFT bin for this target frequency
            % FFT bin k corresponds to frequency k * fs / N
            % So for target frequency f, bin = f * N / fs
            exact_bin = target_freq * data_length / freq;
            
            % Round to nearest integer bin and convert to MATLAB index (+1)
            closest_idx = round(exact_bin) + 1; % +1 for MATLAB 1-based indexing
            
            % Make sure index is within bounds and within positive frequency range
            if closest_idx >= 1 && closest_idx <= nyquist_bin
                harmonics_reindexed(i, col) = harmonics_amplitude(closest_idx, col);
            else
                harmonics_reindexed(i, col) = 0;
            end
        end
    end
    
    % Normalize: divide by length and multiply by 2
    % But DC component (frequency = 0) should only be divided by length
    harmonics_normalized = harmonics_reindexed / data_length * 2;
    
    % Correct DC component normalization (first element if hz_range starts at 0)
    if hz_range(1) == 0
        harmonics_normalized(1, :) = harmonics_reindexed(1, :) / data_length;
    end
    
    % Create output structure
    harmonics = struct();
    harmonics.amplitude = harmonics_normalized;
    harmonics.harmonic = hz_range(:);
    harmonics.type = signal_type;

end
