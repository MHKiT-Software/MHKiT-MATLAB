function result = version_within(actual_version, min_version, max_version, logger)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check if Python version is within specified bounds
%   Checks if actual_version falls within min_version and max_version
%   Only verifies major and minor version bounds
%
% Parameters
% ------------
%     actual_version: string
%         Current Python version (e.g., '3.11.2')
%
%     min_version: string
%         Minimum required version (e.g., '3.10')
%
%     max_version: string
%         Maximum allowed version (e.g., '3.12')
%
%     logger: struct
%         Logger object for error reporting
%
% Returns
% ---------
%     result: logical
%         true if version is within bounds, false otherwise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse all version strings
    actual_parsed = mhkit.python.parse_version_string(actual_version);
    min_parsed = mhkit.python.parse_version_string(min_version);
    max_parsed = mhkit.python.parse_version_string(max_version);

    % Initialize result
    result = true;

    % Check major versions match
    if actual_parsed.major ~= min_parsed.major || actual_parsed.major ~= max_parsed.major
        logger.error('Major versions of python don''t match');
        result = false;
        return;
    end

    % Check minor version bounds
    if actual_parsed.minor < min_parsed.minor || actual_parsed.minor > max_parsed.minor
        logger.error('Python minor version %d is not within bounds [%d, %d]', ...
            actual_parsed.minor, min_parsed.minor, max_parsed.minor);
        result = false;
        return;
    end
end
