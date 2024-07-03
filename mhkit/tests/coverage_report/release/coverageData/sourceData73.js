var sourceData73 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB/mhkit/river/IO/delft_3d/delft_3d_cleanup_turbulent_kinetic_energy.m","RawFileContents":["function result = delft_3d_cleanup_turbulent_kinetic_energy(input, threshold)","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%","% Cleans up the turbulent kinetic energy values based on a threshold.","%","% Parameters","% ------------","%    input: numeric array","%       Array of turbulent kinetic energy values.","%    threshold: numeric","%       Threshold value for cleaning up the turbulent kinetic energy values.","%","% Returns","% ---------","%    result: numeric","%        The cleaned-up turbulent kinetic energy values.","%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","","    if ~isnumeric(input)","        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidInput', 'Input must be numeric.');","    end","","    if ~isnumeric(threshold)","        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidThreshold', 'Threshold must be numeric.');","    end","","    if threshold >= 0","        error('MATLAB:delft_3d_cleanup_turbulent_kinetic_energy:InvalidThresholdValue', 'Threshold must be less than zero.');","    end","","    python_result = py.mhkit_python_utils.delft_3d_helper.cleanup_turbulent_kinetic_energy(py.numpy.array(input), pyargs('threshold', threshold));","    result = double(python_result);","","end",""],"CoverageDisplayDataPerLine":{"Function":{"LineNumber":1,"Hits":1,"StartColumnNumbers":0,"EndColumnNumbers":77,"ContinuedLine":false},"Statement":[{"LineNumber":21,"Hits":1,"StartColumnNumbers":4,"EndColumnNumbers":24,"ContinuedLine":false},{"LineNumber":22,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":105,"ContinuedLine":false},{"LineNumber":25,"Hits":1,"StartColumnNumbers":4,"EndColumnNumbers":28,"ContinuedLine":false},{"LineNumber":26,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":113,"ContinuedLine":false},{"LineNumber":29,"Hits":1,"StartColumnNumbers":4,"EndColumnNumbers":21,"ContinuedLine":false},{"LineNumber":30,"Hits":0,"StartColumnNumbers":8,"EndColumnNumbers":125,"ContinuedLine":false},{"LineNumber":33,"Hits":1,"StartColumnNumbers":4,"EndColumnNumbers":146,"ContinuedLine":false},{"LineNumber":34,"Hits":1,"StartColumnNumbers":4,"EndColumnNumbers":35,"ContinuedLine":false}]}}