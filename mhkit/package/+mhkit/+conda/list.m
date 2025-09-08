function packages = list(env_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gets a list of packages from conda environment
%
% Parameters
% ------------
%     env_name: string
%         Conda environment name
%
% Returns
% ---------
%     packages: struct
%         Structure where each field is a package name containing a struct
%         with 'version', 'build', and 'channel' fields
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   disp(packages.numpy.version);
%
%   % Parse packages from a file
%   packages = get_conda_packages('conda_list.txt', true);
%   disp(packages.pip.version);

% Default is_file to false if not provided
if nargin < 2
    is_file = false;
end

% Initialize the output struct
packages = struct();

% Get conda output either from command or file
% Run conda command
% Construct command to activate environment and list packages
cmd = sprintf('conda run -n %s conda list', env_name);

[status, cmdout] = mhkit.sys(cmd);

% Check if command was successful
if status ~= 0
    error('Failed to execute conda command. Error: %s', cmdout);
end

conda_output = cmdout;

% Split into lines
lines = strsplit(conda_output, '\n');

% Find the line that starts with '#' followed by spaces and 'Name'
header_line_idx = find(contains(lines, '# Name'), 1);

% If no header found, return empty struct
if isempty(header_line_idx)
    warning('Could not find package header in conda output');
    return;
end

% Process each line after the header line
for i = (header_line_idx+2):length(lines) % +2 to skip header and the line with dashes
    line = strtrim(lines{i});

    % Skip empty lines and comments
    if isempty(line) || startsWith(line, '#')
        continue;
    end

    % Split the line by whitespace
    parts = regexp(line, '\s+', 'split');

    % Skip if we don't have at least name and version
    if length(parts) < 2
        continue;
    end

    % Get package name and make it a valid field name
    pkg_name = parts{1};

    % Some package names might contain characters that aren't valid MATLAB field names
    % Replace invalid characters with underscores
    valid_field_name = matlab.lang.makeValidName(pkg_name);

    % Create package info struct
    pkg_info = struct();
    pkg_info.version = parts{2};

    % Store original name in case it was modified
    pkg_info.name = pkg_name;

    % Add build and channel if available
    if length(parts) >= 3
        pkg_info.build = parts{3};
    else
        pkg_info.build = '';
    end

    if length(parts) >= 4
        pkg_info.channel = parts{4};
    else
        pkg_info.channel = '';
    end

    % Add to the output struct
    packages.(valid_field_name) = pkg_info;
end

% If no packages were found, display a warning
if isempty(fieldnames(packages))
    warning('No packages found in the conda output');
end

end
