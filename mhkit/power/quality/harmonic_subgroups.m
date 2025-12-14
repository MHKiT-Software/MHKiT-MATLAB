function harmonic_subgroups = harmonic_subgroups(harmonics, grid_freq)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the harmonic subgroups based on IEC 61000-4-7
%     
% Computes harmonic subgroups by grouping adjacent harmonics around 
% target frequencies determined by the grid frequency.
%
% Parameters
% ------------
%   harmonics: structure
%       harmonics.amplitude : Harmonic amplitude data with each timeseries in its own column
%       harmonics.harmonic : Harmonic frequency vector [Hz]
%   grid_freq: numeric
%       Value indicating if the power supply is 50 or 60 Hz. Valid inputs are 50 and 60
%
% Returns
% ---------
%   harmonic_subgroups: structure
%       harmonic_subgroups.amplitude : Harmonic subgroup amplitudes
%       harmonic_subgroups.harmonic : Target harmonic frequencies [Hz]
%
% Key Equations
% -------------
% 1. Subgroup calculation:
%    subgroup_value = sqrt(h_prev^2 + h_current^2 + h_next^2)
%    where h represents harmonic amplitudes at adjacent frequencies
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStruct = @(x) isstruct(x);
    validGridFreq = @(x) isnumeric(x) && (x == 50 || x == 60);
    
    % Add required parameters
    addRequired(p, 'harmonics', validStruct);
    addRequired(p, 'grid_freq', validGridFreq);
    
    % Parse inputs
    parse(p, harmonics, grid_freq);
    
    % Extract validated inputs
    harmonics = p.Results.harmonics;
    grid_freq = p.Results.grid_freq;
    
    % Validate input structures have required fields
    if ~isfield(harmonics, 'amplitude')
        error('MHKiT:harmonic_subgroups: harmonics structure must contain amplitude field');
    end
    if ~isfield(harmonics, 'harmonic')
        error('MHKiT:harmonic_subgroups: harmonics structure must contain harmonic field');
    end
    
    % Extract data
    harmonic_data = harmonics.amplitude;
    harmonic_freq = harmonics.harmonic;
    
    % Validate dimensions
    if size(harmonic_data, 1) ~= length(harmonic_freq)
        error('MHKiT:harmonic_subgroups: Number of harmonic frequencies must match the number of rows in amplitude data');
    end
    
    % Validate grid frequency
    if ~(grid_freq == 50 || grid_freq == 60)
        error('MHKiT:harmonic_subgroups: grid_freq must be either 50 or 60');
    end
    
    % Create target frequency array based on grid frequency
    if grid_freq == 60
        target_hz = (0:60:3000)';
    else % grid_freq == 50
        target_hz = (0:50:2500)';
    end
    
    % Sort harmonics by frequency for proper indexing
    [sorted_freq, sort_idx] = sort(harmonic_freq);
    sorted_data = harmonic_data(sort_idx, :);
    
    % Get number of columns (signals)
    num_signals = size(harmonic_data, 2);
    num_targets = length(target_hz);
    
    % Initialize output array
    subgroup_amplitudes = zeros(num_targets, num_signals);
    
    % Calculate harmonic subgroups
    for i = 1:num_targets
        current_target = target_hz(i);
        
        % Find nearest harmonic frequency index
        [~, nearest_idx] = min(abs(sorted_freq - current_target));
        
        % Define indices for the three adjacent harmonics (with bounds checking)
        idx_prev = max(1, nearest_idx - 1);
        idx_curr = nearest_idx;
        idx_next = min(length(sorted_freq), nearest_idx + 1);
        
        % Calculate RMS sum for each signal column
        for j = 1:num_signals
            prev_val = sorted_data(idx_prev, j);
            curr_val = sorted_data(idx_curr, j);
            next_val = sorted_data(idx_next, j);
            
            % Calculate subgroup as RMS sum of three adjacent harmonics
            subgroup_amplitudes(i, j) = sqrt(prev_val^2 + curr_val^2 + next_val^2);
        end
    end
    
    % Create output structure
    harmonic_subgroups = struct();
    harmonic_subgroups.amplitude = subgroup_amplitudes;
    harmonic_subgroups.harmonic = target_hz;

end
