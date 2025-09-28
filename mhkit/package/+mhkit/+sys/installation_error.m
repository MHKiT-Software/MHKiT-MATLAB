function installation_error(logger, error_message, step_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Standard error handler for installation failures
%
% Parameters
% ------------
%     logger: struct
%         Logger object for outputting status messages
%     error_message: string
%         Specific error message describing what failed
%     step_name: string
%         Name of the installation step that failed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spec = mhkit.spec();

    logger.error('%s', error_message);
    logger.error('');
    logger.error('Installation step "%s" failed.', step_name);
    logger.error('');
    logger.error('Next steps:');
    logger.error('  1. Check system requirements and network connection');
    if ispc
        logger.error('  2. For conda issues: %s', spec.support.windows_conda_install);
    else
        logger.error('  2. For conda issues: %s', spec.support.conda_install_instructions);
    end
    logger.error('  3. Report this issue: %s', spec.support.github_issues);
    logger.error('  4. Include error details and system information when reporting');
    logger.error('');
end