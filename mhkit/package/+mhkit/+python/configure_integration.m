function configure_integration(python_path, execution_mode, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure MATLAB's Python integration using pyenv
%
% Parameters
% ------------
%     python_path: string
%         Full path to the Python executable
%     execution_mode: string
%         MATLAB Python execution mode ('InProcess' or 'OutOfProcess')
%     logger: struct
%         Logger object for outputting status messages
%
% Returns
% ---------
%     None (throws error on failure)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    logger.info('Configuring MATLAB Python integration...');
    logger.info('  Python path: %s', python_path);
    logger.info('  Execution mode: %s', execution_mode);

    % Add Python directory to system PATH
    python_dir = fileparts(python_path);
    current_path = getenv('PATH');
    if isunix
        new_path = [python_dir ':' current_path];
    else
        new_path = [python_dir ';' current_path];
    end
    setenv('PATH', new_path);
    logger.info('✓ Added Python directory to PATH: %s', python_dir);

    % Configure pyenv
    try
        logger.info('Configuring pyenv...');
        pyenv(Version=python_path, ExecutionMode=execution_mode);
        logger.info('✓ pyenv command executed successfully');

    catch ME
        error_msg = sprintf('pyenv configuration failed. Python path: %s, ExecutionMode: %s. Error: %s', ...
                          python_path, execution_mode, ME.message);
        mhkit.sys.installation_error(logger, error_msg, 'Python Integration');
        error(error_msg);
    end

    % Verify pyenv configuration
    logger.info('Verifying pyenv configuration...');
    try
        pe = pyenv();
        logger.info('pyenv status:');
        logger.info('  Version: "%s"', pe.Version);
        logger.info('  Executable: "%s"', pe.Executable);
        logger.info('  Status: "%s"', pe.Status);
        logger.info('  ExecutionMode: "%s"', pe.ExecutionMode);

        % Check if execution mode was set correctly
        if ~strcmp(pe.ExecutionMode, execution_mode)
            logger.warning('ExecutionMode mismatch: expected "%s", got "%s"', execution_mode, pe.ExecutionMode);
        end

        % Check if Python is configured
        if strcmp(pe.Status, 'NotLoaded')
            logger.info('Python status is "NotLoaded" - this is normal before first Python call');
        else
            logger.info('✓ Python status: %s', pe.Status);
        end

    catch ME
        error_msg = sprintf('Failed to verify pyenv configuration: %s', ME.message);
        mhkit.sys.installation_error(logger, error_msg, 'Python Integration');
        error(error_msg);
    end

    logger.info('✓ Python integration configured successfully');
end