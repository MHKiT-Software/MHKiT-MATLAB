function analytic_signal = mhkit_scipy_hilbert(signal)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates Hilbert transform using SciPy's implementation
%
% Parameters
% ------------
%   signal: array
%       Input signal for Hilbert transform
%
% Returns
% ---------
%   analytic_signal: array (complex)
%       Analytic signal from Hilbert transform
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if SciPy is available
    if ~mhkit_check_scipy()
        error('MHKiT: scipy not available. Please follow the documentation instructions to install python and scipy: https://mhkit-software.github.io/MHKiT/matlab_installation.html');
    end
    
    try
        % Convert MATLAB array to Python
        py_signal = py.numpy.array(signal);
        
        % Call SciPy's Hilbert transform
        py_result = py.scipy.signal.hilbert(py_signal);
        
        % Convert back to MATLAB
        analytic_signal = double(py_result);
        
    catch ME
        error('MHKiT: Error calling SciPy hilbert: %s. Please report this issue: https://github.com/MHKiT-Software/MHKiT-MATLAB/issues', ME.message);
    end
end
