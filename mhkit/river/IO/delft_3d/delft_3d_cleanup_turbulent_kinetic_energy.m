function result = delft_3d_cleanup_turbulent_kinetic_energy(input, threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cleans up the turbulent kinetic energy values based on a threshold.
%
% Parameters
% ------------
%    input: numeric array
%       Array of turbulent kinetic energy values.
%    threshold: numeric
%       Threshold value for cleaning up the turbulent kinetic energy values.
%
% Returns
% ---------
%    result: numeric
%        The cleaned-up turbulent kinetic energy values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isnumeric(input)
        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidInput', 'Input must be numeric.');
    end

    if ~isnumeric(threshold)
        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidThreshold', 'Threshold must be numeric.');
    end

    if threshold >= 0
        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidThresholdValue', 'Threshold must be less than zero.');
    end

    python_result = py.mhkit_python_utils.delft_3d_helper.cleanup_turbulent_kinetic_energy(py.numpy.array(input), pyargs('threshold', threshold));
    result = double(python_result);

end
