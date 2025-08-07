function is_installed = mhkit_check_scipy()
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Checks if Python is installed and SciPy is accessible from MATLAB
%
% Returns
% ---------
%   is_installed: logical
%       true if both Python and SciPy are available, false otherwise
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    is_installed = false;  % Default to false
    
    % First check if Python is available
    try
        % Try to get Python version
        version_info = pyversion;
        if isempty(version_info) || strcmp(version_info, '')
            % Silent return - Python not available
            return;
        end
    catch
        % Silent return - Python check failed
        return;
    end
    
    % Now check if SciPy is available
    try
        % Try to import scipy
        scipy_module = py.importlib.import_module('scipy');
        
        % Try to actually use scipy.signal (what we need)
        test_signal = py.numpy.array([1, 2, 3, 4, 5]);
        py.scipy.signal.hilbert(test_signal);
        
        % If we get here, both scipy and scipy.signal.hilbert work
        is_installed = true;
        
    catch
        % Silent return - SciPy not available or hilbert function failed
        return;
    end
end
