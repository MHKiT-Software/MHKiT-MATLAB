function THCD = total_harmonic_current_distortion(harmonic_subgroups, rated_current)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the total harmonic current distortion (THC) based on IEC TS 62600-30
%
% Parameters
% ------------
%   harmonic_subgroups: structure
%       harmonic_subgroups.amplitude : Subgrouped current harmonics amplitude indexed by harmonic order [Amps]
%       harmonic_subgroups.harmonic : Harmonic frequency order vector
%   rated_current: double
%       Rated current of the energy device [Amps]
%
% Returns
% ---------
%   THCD: double
%       Total harmonic current distortion [%]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments (Input)
        harmonic_subgroups struct
        rated_current (1,1) {mustBeNumeric, mustBePositive}
    end

    arguments (Output)
        THCD
    end
    
    % Validate input structures have required fields
    if ~isfield(harmonic_subgroups, 'amplitude')
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'harmonic_subgroups structure must contain amplitude field');
    end
    if ~isfield(harmonic_subgroups, 'harmonic')
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'harmonic_subgroups structure must contain harmonic field');
    end
    
    % Extract amplitude data
    harmonics_data = harmonic_subgroups.amplitude;
    
    % Validate amplitude data dimensions
    if ~isnumeric(harmonics_data)
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'harmonic_subgroups.amplitude must be numeric');
    end
    
    % Check if we have enough harmonic data (need at least index 2 for fundamental and some harmonics)
    if length(harmonics_data) < 3
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'harmonic_subgroups.amplitude must contain at least 3 elements');
    end
    
    % Convert Python indexing to MATLAB indexing
    % Python: harmonics_subgroup.iloc[2:50] means indices 2 through 49 (0-based)
    % MATLAB: equivalent is indices 3 through 50 (1-based)
    harmonic_start_idx = 3;  % Python index 2 + 1
    harmonic_end_idx = min(50, length(harmonics_data));  % Python index 49 + 1, but limit to data length
    
    % Extract harmonic components (excluding fundamental)
    % Python: harmonics_subgroup.iloc[2:50]**2
    % MATLAB: harmonics_data(3:end_idx).^2 (element-wise power)
    if harmonic_end_idx >= harmonic_start_idx
        harmonics_subset = harmonics_data(harmonic_start_idx:harmonic_end_idx);
    else
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'Insufficient harmonic data for calculation');
    end
    
    % Square the harmonic amplitudes (element-wise operation)
    harmonics_sq = harmonics_subset .^ 2;
    
    % Sum the squared harmonics
    harmonics_sum = sum(harmonics_sq);
    
    % Get fundamental component
    % Python: harmonics_subgroup.iloc[1] means index 1 (0-based)
    % MATLAB: harmonics_data(2) means index 2 (1-based)
    if length(harmonics_data) < 2
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'harmonic_subgroups.amplitude must contain fundamental component at index 2');
    end
    
    fundamental_current = harmonics_data(2);  % Python index 1 + 1
    
    % Validate fundamental current is not zero
    if fundamental_current == 0
        error('MHKiT:total_harmonic_current_distortion:InvalidInput', 'Fundamental current component cannot be zero');
    end
    
    % Calculate Total Harmonic Current Distortion
    % Python: (np.sqrt(harmonics_sum)/harmonics_subgroup.iloc[1])*100
    % MATLAB: (sqrt(harmonics_sum) ./ fundamental_current) .* 100
    THCD = (sqrt(harmonics_sum) ./ fundamental_current) .* 100;
    
    % Ensure output is a scalar double
    THCD = double(THCD);

end
