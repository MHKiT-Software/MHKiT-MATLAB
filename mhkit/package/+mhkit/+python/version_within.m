function result = version_within(actual_version, min_version, max_version, logger)
    % VERSION_WITHIN Check if Python version is within specified bounds
    %   Checks if actual_version falls within min_version and max_version
    %   Only verifies major and minor version bounds
    %
    %   Input:
    %       actual_version - string, current Python version
    %       min_version - string, minimum required version
    %       max_version - string, maximum allowed version
    %       logger - logger object for error reporting
    %
    %   Output:
    %       result - boolean, true if version is within bounds

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
