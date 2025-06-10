function interharmonics = interharmonics(harmonics, grid_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the interharmonics from the harmonics of current
%     
% This function computes interharmonic groups by analyzing harmonic amplitude data
% at specific frequency intervals.
%
% Parameters
% ------------
%   harmonics: structure
%       harmonics.amplitude : Harmonic amplitude data with each timeseries in its own column [dimensionless]
%       harmonics.harmonic : Harmonic frequency vector corresponding to amplitude data [Hz]
%   grid_freq: numeric scalar
%       Value indicating if the power supply is 50 or 60 Hz. Options = 50 or 60 [Hz]
%
% Returns
% ---------
%   interharmonics: structure
%       interharmonics.amplitude : Interharmonic group amplitudes [same units as input]
%       interharmonics.harmonic : Frequency vector for interharmonic groups [Hz]
%
% Key Equations
% -------------
% 1. Interharmonic group calculation:
%    IH_group = sqrt(sum(H_subset^2))
%    where H_subset contains harmonics between consecutive grid frequencies
%
% 2. Frequency grid generation:
%    For 60 Hz: f_grid = 0:60:3060
%    For 50 Hz: f_grid = 0:50:2550
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStruct = @(x) isstruct(x);
    validGridFreq = @(x) isnumeric(x) && isscalar(x) && (x == 50 || x == 60);
    
    % Add required parameters
    addRequired(p, 'harmonics', validStruct);
    addRequired(p, 'grid_freq', validGridFreq);
    
    % Parse inputs
    parse(p, harmonics, grid_freq);
    
    % Extract validated inputs
    harmonics = p.Results.harmonics;
    grid_freq = p.Results.grid_freq;
    
    % Validate input structure has required fields
    if ~isfield(harmonics, 'amplitude')
        error('MHKiT:interharmonics: harmonics structure must contain amplitude field');
    end
    if ~isfield(harmonics, 'harmonic')
        error('MHKiT:interharmonics: harmonics structure must contain harmonic field');
    end
    
    % Extract data from structure
    harmonics_amplitude = harmonics.amplitude;
    harmonics_frequency = harmonics.harmonic;
    
    % Validate grid frequency
    if ~(grid_freq == 50 || grid_freq == 60)
        error('MHKiT:interharmonics: grid_freq must be either 50 or 60');
    end
    
    % Validate dimensions
    if size(harmonics_amplitude, 1) ~= length(harmonics_frequency)
        error('MHKiT:interharmonics: harmonics.amplitude rows must match length of harmonics.harmonic');
    end
    
    % Create frequency grid based on grid frequency
    if grid_freq == 60
        hz_grid = 0:60:3060;  % MATLAB: 1-based, but this creates the correct range
        hz_grid = hz_grid(1:end-1);  % Remove last element to match Python behavior
        subset_size = 10;  % Number of harmonics to include in each interharmonic group
    elseif grid_freq == 50
        hz_grid = 0:50:2550;  % MATLAB: 1-based, but this creates the correct range  
        hz_grid = hz_grid(1:end-1);  % Remove last element to match Python behavior
        subset_size = 6;   % Number of harmonics to include in each interharmonic group
    end
    
    % Get number of columns in amplitude data
    num_cols = size(harmonics_amplitude, 2);
    num_frequencies = length(hz_grid);
    
    % Initialize output matrix
    interharmonics_amplitude = ones(num_frequencies, num_cols);
    
    % Sort harmonics by frequency (equivalent to sort_index in Python)
    [sorted_freq, sort_idx] = sort(harmonics_frequency);
    sorted_amplitude = harmonics_amplitude(sort_idx, :);
    
    % Process each frequency in the grid
    for i = 1:num_frequencies  % MATLAB: 1-based indexing
        current_freq = hz_grid(i);
        
        % Find nearest frequency index (equivalent to get_loc with method='nearest')
        [~, nearest_idx] = min(abs(sorted_freq - current_freq));
        
        % Process each column of amplitude data
        for j = 1:num_cols  % MATLAB: 1-based indexing
            
            % Calculate subset indices with bounds checking
            % Python: harmonics[col].iloc[indn+1:indn+11] (60Hz) or [indn+1:indn+7] (50Hz)
            % MATLAB: Need to add 1 for 1-based indexing and handle bounds
            start_idx = nearest_idx + 1;  % Python indn+1 becomes nearest_idx+1 in MATLAB
            end_idx = start_idx + subset_size - 1;  % Python slice end is exclusive, MATLAB inclusive
            
            % Ensure indices are within bounds
            if start_idx > length(sorted_freq)
                % If start is beyond array, use zeros
                subset_squared = zeros(subset_size, 1);
            elseif end_idx > length(sorted_freq)
                % If end is beyond array, take what we can and pad with zeros
                available_size = length(sorted_freq) - start_idx + 1;
                if available_size > 0
                    subset_data = sorted_amplitude(start_idx:end, j);
                    subset_squared = [subset_data.^2; zeros(subset_size - available_size, 1)];
                else
                    subset_squared = zeros(subset_size, 1);
                end
            else
                % Normal case: extract subset and square
                subset_data = sorted_amplitude(start_idx:end_idx, j);
                subset_squared = subset_data.^2;  % MATLAB: .^ for element-wise power
            end
            
            % Calculate interharmonic as square root of sum of squares
            interharmonics_amplitude(i, j) = sqrt(sum(subset_squared));
        end
    end
    
    % Create output structure
    interharmonics = struct();
    interharmonics.amplitude = interharmonics_amplitude;
    interharmonics.harmonic = hz_grid';  % Convert to column vector to match input format

end
