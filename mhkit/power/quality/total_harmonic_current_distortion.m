function THCD = total_harmonic_current_distortion(harmonic_subgroups, rated_current)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the total harmonic current distortion (THC) based on IEC/TS 62600-30
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
% Key Equations
% -------------
% 1. Harmonic current distortion calculation:
%    THCD = (sqrt(sum(I_h^2)) / I_1) * 100
%    where I_h are harmonic currents (orders 2-50) and I_1 is fundamental current
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create input parser
    p = inputParser;
    
    % Define validation functions
    validStruct = @(x) isstruct(x);
    validNumeric = @(x) isnumeric(x) && isscalar(x);
    
    % Add required parameters
    addRequired(p, 'harmonic_subgroups', validStruct);
    addRequired(p, 'rated_current', validNumeric);
    
    % Parse inputs
    parse(p, harmonic_subgroups, rated_current);
    
    % Extract validated inputs
    harmonic_subgroups = p.Results.harmonic_subgroups;
    rated_current = p.Results.rated_current;
    
    % Validate input structures have required fields
    if ~isfield(harmonic_subgroups, 'amplitude')
        error('MHKiT:total_harmonic_current_distortion: harmonic_subgroups structure must contain amplitude field');
    end
    if ~isfield(harmonic_subgroups, 'harmonic')
        error('MHKiT:total_harmonic_current_distortion: harmonic_subgroups structure must contain harmonic field');
    end
    
    % Extract amplitude data
    harmonics_data = harmonic_subgroups.amplitude;
    
    % Validate amplitude data dimensions
    if ~isnumeric(harmonics_data)
        error('MHKiT:total_harmonic_current_distortion: harmonic_subgroups.amplitude must be numeric');
    end
    
    % Check if we have enough harmonic data (need at least index 2 for fundamental and some harmonics)
    if length(harmonics_data) < 3
        error('MHKiT:total_harmonic_current_distortion: harmonic_subgroups.amplitude must contain at least 3 elements');
    end
    
    % Validate rated_current is positive
    if rated_current <= 0
        error('MHKiT:total_harmonic_current_distortion: rated_current must be positive');
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
        error('MHKiT:total_harmonic_current_distortion: Insufficient harmonic data for calculation');
    end
    
    % Square the harmonic amplitudes (element-wise operation)
    harmonics_sq = harmonics_subset .^ 2;
    
    % Sum the squared harmonics
    harmonics_sum = sum(harmonics_sq);
    
    % Get fundamental component
    % Python: harmonics_subgroup.iloc[1] means index 1 (0-based)
    % MATLAB: harmonics_data(2) means index 2 (1-based)
    if length(harmonics_data) < 2
        error('MHKiT:total_harmonic_current_distortion: harmonic_subgroups.amplitude must contain fundamental component at index 2');
    end
    
    fundamental_current = harmonics_data(2);  % Python index 1 + 1
    
    % Validate fundamental current is not zero
    if fundamental_current == 0
        error('MHKiT:total_harmonic_current_distortion: Fundamental current component cannot be zero');
    end
    
    % Calculate Total Harmonic Current Distortion
    % Python: (np.sqrt(harmonics_sum)/harmonics_subgroup.iloc[1])*100
    % MATLAB: (sqrt(harmonics_sum) ./ fundamental_current) .* 100
    THCD = (sqrt(harmonics_sum) ./ fundamental_current) .* 100;
    
    % Ensure output is a scalar double
    THCD = double(THCD);

end
