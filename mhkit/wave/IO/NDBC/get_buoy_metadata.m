function data = get_buoy_metadata(station_number)
% GET_BUOY_METADATA Fetches and parses metadata for a National Data Buoy Center (NDBC) station.
%
% This function extracts information such as provider, buoy type, latitude, longitude,
% and other metadata from the station's webpage.
%
% Parameters
% ----------
% station_number : string
%     The station number (ID) of the NDBC buoy.
%
% Returns
% -------
% data : struct
%     A structure containing metadata of the buoy with fields representing
%     different types of information.
%
% Example
% -------
%     metadata = get_buoy_metadata('44025');
%     disp(metadata);
%
% Dependencies
% ------------
% - MATLAB's built-in web functions (`webread`, `strfind`, etc.)
% - Requires the HTML source of the station page to be available.

% Define the URL for the station
url = sprintf('https://www.ndbc.noaa.gov/station_page.php?station=%s', station_number);

try
    % Fetch the page content
    content = webread(url);
catch
    error('Failed to fetch data. Check the station number or network connection.');
end

% Initialize output data structure
data = struct();

% Check if the station is valid
if contains(content, 'Station not found')
    error('Invalid or nonexistent station number: %s', station_number);
end

% Extract the station title (buoy name)
title_start = strfind(content, '<h1>') + 4;
title_end = strfind(content, '</h1>') - 1;
if isempty(title_start) || isempty(title_end)
    error('Failed to extract the buoy title.');
end
title = strtrim(content(title_start:title_end));
title = regexprep(title, '\n.*$', '');  % Remove any trailing newlines and extra content
data.buoy = title;

% Find the specific div containing the buoy metadata
metadata_start = strfind(content, '<div id="stn_metadata">');
if isempty(metadata_start)
    warning('Metadata section not found.');
    return;
end
metadata_end = strfind(content(metadata_start:end), '</div>') + metadata_start - 1;
metadata_text = content(metadata_start:metadata_end);

% Split metadata into lines
lines = splitlines(metadata_text);

% Parse each line
line_count = 1;
for i = 1:length(lines)
    line = strtrim(lines{i});

    % Check for latitude/longitude pattern
    if ~isempty(regexp(line, '\d+\.\d+\s+[NS]\s+\d+\.\d+\s+[EW]', 'once'))
        tokens = regexp(line, '(\d+\.\d+)\s+([NS])\s+(\d+\.\d+)\s+([EW])', 'tokens');
        lat = tokens{1}{1};
        lon = tokens{1}{3};
        data.lat = lat;
        data.lon = lon;
    elseif contains(line, ':')
        % Split key-value pairs on colon
        parts = strsplit(line, ':');
        key = strtrim(parts{1});
        value = strtrim(parts{2});
        data.(matlab.lang.makeValidName(key)) = value;
    elseif ~isempty(line)
        % Catch all other lines as keys with empty values
        data.(matlab.lang.makeValidName(line)) = '';
    end
    line_count = line_count + 1;
end
end
