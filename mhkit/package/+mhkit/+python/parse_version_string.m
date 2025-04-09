function version = parse_python_version(version_str)
    % PARSE_PYTHON_VERSION Parse Python version string into a structured format
    %   version = PARSE_PYTHON_VERSION(version_str) takes a Python version string
    %   and returns a struct with major, minor, patch, and additional version info
    %
    %   Input:
    %       version_str - string, Python version (e.g., "3.11.8.final.0")
    %
    %   Output:
    %       version - struct containing:
    %           - major (integer)
    %           - minor (integer)
    %           - patch (integer)
    %           - qualifier (string, e.g., 'final')
    %           - build (integer)
    %           - original (string, the original version string)
    %
    %   Example:
    %       ver = parse_python_version("3.11.8.final.0")
    %       % Returns:
    %       %   ver.major = 3
    %       %   ver.minor = 11
    %       %   ver.patch = 8
    %       %   ver.qualifier = 'final'
    %       %   ver.build = 0
    %       %   ver.original = "3.11.8.final.0"

    % Initialize the output struct
    disp("Parsing python version from:");
    disp(version_str);
    version = struct();
    version.original = version_str;

    % Split the version string
    parts = strsplit(version_str, '.');

    % Extract major.minor.patch
    version.major = str2double(parts{1});
    version.minor = str2double(parts{2});

    if length(parts) > 2
        version.patch = str2double(parts{3});
    else
        version.patch = 0
    end

    % Extract qualifier and build if they exist
    if length(parts) > 3
        version.qualifier = parts{4};
    else
        version.qualifier = '';
    end

    if length(parts) > 4
        version.build = str2double(parts{5});
    else
        version.build = 0;
    end

    % Validate numeric fields
    if isnan(version.major) || isnan(version.minor) || isnan(version.patch)
        error('Invalid version format: major, minor, and patch must be numeric');
    end

    % Ensure all numeric fields are integers
    version.major = floor(version.major);
    version.minor = floor(version.minor);
    version.patch = floor(version.patch);
    version.build = floor(version.build);
end
