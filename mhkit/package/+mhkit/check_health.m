function result = check_health()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Validate MHKiT Python integration and environment health
%   Performs comprehensive checks of Python environment, module imports,
%   and MHKiT functionality to ensure everything is working correctly
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     result: logical
%         true if all health checks pass, false if any fail
%
% Example
% -------
%     if mhkit.check_health()
%         disp('MHKiT is ready to use!');
%     else
%         disp('MHKiT has issues - check the output above');
%     end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    logger = mhkit.logging.get_logger();
    result = true; % Start optimistic
    
    try
        logger.info('=== MHKiT Health Check ===');
        
        % Check 1: Python Environment Configuration
        logger.info('=== Python Environment Validation ===');
        try
            pe = pyenv();
            logger.info('Python Version: "%s"', pe.Version);
            logger.info('Python Executable: "%s"', pe.Executable);
            logger.info('Python Library: "%s"', pe.Library);
            logger.info('Python Home: "%s"', pe.Home);
            logger.info('Execution Mode: "%s"', pe.ExecutionMode);
            logger.info('Status: "%s"', pe.Status);
            
            % Additional debugging for empty fields
            if isempty(pe.Version)
                logger.warning('  Python Version is empty - pyenv not configured');
            end
            if isempty(pe.Executable)
                logger.warning('  Python Executable is empty - no Python path set');
            end
            
            % Validate Python is properly configured
            if strcmp(pe.Status, 'NotLoaded')
                logger.error('✗ Python is not loaded');
                result = false;
            elseif strcmp(pe.Status, 'Loaded')
                logger.info('✓ Python environment is loaded');
            else
                logger.warning('  Python status is: %s', pe.Status);
            end
            
        catch ME
            logger.error('✗ Failed to get Python environment info: %s', ME.message);
            result = false;
        end
        
        % Check 2: Python Module Import Tests
        logger.info('=== Testing Python Module Import ===');
        
        % Test mhkit module import
        try
            py.importlib.import_module('mhkit');
            logger.info('✓ mhkit module imported successfully');
        catch ME
            logger.error('✗ Failed to import mhkit: %s', ME.message);
            logger.error('  Try: pip install mhkit==0.9');
            result = false;
        end
        
        % Test mhkit_python_utils module import
        try
            py.importlib.import_module('mhkit_python_utils');
            logger.info('✓ mhkit_python_utils module imported successfully');
        catch ME
            logger.error('✗ Failed to import mhkit_python_utils: %s', ME.message);
            logger.error('  This module should be installed by mhkit.install()');
            result = false;
        end
        
        % Check 3: MHKiT Functionality Test
        logger.info('=== Testing MHKiT Functionality ===');
        try
            test_result = py.mhkit.river.performance.circular(30);
            logger.info('✓ MHKiT function test successful: %s', char(test_result));
        catch ME
            logger.error('✗ Failed to execute MHKiT function: %s', ME.message);
            logger.error('  This indicates a problem with the MHKiT Python installation');
            result = false;
        end
        
        % Check 4: Additional System Information
        logger.info('=== System Information ===');
        try
            % Check MATLAB version
            logger.info('MATLAB Version: %s', version('-release'));
            
            % Check platform
            if ispc
                platform = 'Windows';
            elseif ismac
                platform = 'macOS';
            elseif isunix
                platform = 'Linux';
            else
                platform = 'Unknown';
            end
            logger.info('Platform: %s', platform);
            
            % Check PATH for Python
            current_path = getenv('PATH');
            if contains(current_path, 'python') || contains(current_path, 'conda')
                logger.info('✓ Python/Conda found in system PATH');
            else
                logger.warning('⚠ Python/Conda not detected in system PATH');
            end
            
        catch ME
            logger.warning('⚠ Could not gather all system information: %s', ME.message);
        end
        
        % Final Result
        if result
            logger.info('=== ✓ All Health Checks Passed ===');
            logger.info('MHKiT is ready to use!');
        else
            logger.error('=== ✗ Health Check Failed ===');
            logger.error('MHKiT has configuration issues that need to be resolved.');
            logger.error('');
            logger.error('Troubleshooting steps:');
            logger.error('1. Try restarting MATLAB');
            logger.error('2. Run mhkit.install() again');
            logger.error('3. Check that pyenv() points to the correct Python environment');
            logger.error('4. Verify conda/Python are in your system PATH');
        end
        
    catch ME
        logger.error('=== ✗ Health Check Error ===');
        logger.error('Health check encountered an unexpected error: %s', ME.message);
        result = false;
    end
end
