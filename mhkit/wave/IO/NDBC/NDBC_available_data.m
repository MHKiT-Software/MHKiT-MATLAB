function available_data=NDBC_available_data(parameter,options)

%%%%%%%%%%%%%%%%%%%%
%     For a given parameter this will return a structure of years, 
%     station IDs and file names that contain that parameter data. 
%     
%     
% Parameters
% ------------
%     parameter : string
%         'swden' : 'Raw Spectral Wave Current Year Historical Data'
%         'stdmet': 'Standard Meteorological Current Year Historical Data'
%
%     buoy_number : string (optional)
%         Buoy Number.  5-character alpha-numeric station identifier 
%         to call: NDBC_available_data(parameter,"buoy_number",buoy_number)
%
%     proxy : None
%         Parameter is now deprecated.
%         To request data from behind a firewall, configure in MATLAB
%         Preferences by navigating to:
%           Home -> Environment -> Preferences
%         then:
%           MATLAB -> Web -> Use a proxy server to connect to the Internet
%         See the following for details:
%           https://www.mathworks.com/help/matlab/import_export/proxy.html
%     
% Returns
% ---------
%     available_data : structure 
%         Structure with station_id, years, and NDBC file names.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    parameter string {mustBeMember(parameter,{'swden','stdmet'})}
    options.buoy_number string = "";
end

if strlength(options.buoy_number) ~= 0 && strlength(options.buoy_number) ~= 5
    ME = MException('MATLAB:NDBC_available_data', ...
                    'Buoy must be a 5-character alpha-numeric identifier.');
end

MAX_RETRIES = 5;                         % number of query retries if error
FILENAME_IDENTIFIER = ".txt.gz";

% Formulate query
data_url = "https://www.ndbc.noaa.gov/data/historical/";
api_query = parameter;

% Display query
disp("Data request URL: " + data_url + api_query)

% Submit query and get data
for i = 0:MAX_RETRIES
    try
        response = webread(data_url + api_query);
        break;
    catch ME
        if i == MAX_RETRIES
            disp(['MATLAB:NDBC_available_data: ', ME.identifier]);
            rethrow(ME)
        else
            pause(1);   % pause(seconds) and retry query
        end
    end
end

% Organize a structure containing lists of the available data
data_lines = strsplit(response, '\n');
available_data.Station_id = [];
available_data.year = [];
available_data.file = [];
for i = 1:length(data_lines)
    if contains(data_lines{i}, FILENAME_IDENTIFIER)
        % Parse and add filename to output structure
        % Note: need to grab between '"' as escaping '>' or '<' has issues
        % Example filename: '42otpw2000.txt.gz', where:
        %   2000 = year
        %   w = delimiter (can also be an 'h')
        regex_pattern = strcat( ...
            '(?<=\")', ...       % look after a '"'
            '[a-zA-Z0-9_]*', ... % any number of these alphanumerics and underscore
            strrep(FILENAME_IDENTIFIER, '.', '\.'), ...  % escape the periods
            '(?=\")' ...         % look before a '"'
            );
        filename = regexp(data_lines{i}, regex_pattern, 'match');
        available_data.file = [available_data.file; ...
                               convertCharsToStrings(filename{1})];

        % Parse year from filename and add to output structure
        filename_no_extension = extractBefore(filename, '.');
        if contains(filename_no_extension, '_')
            % Handles filenames like: 42002w2008_old.txt.gz
            filename_base = extractBefore(filename_no_extension, '_');
        else
            filename_base = filename_no_extension;
        end
        year = str2double(extractAfter( ...
            filename_base, ...
            strlength(filename_base) - 4));     % grab last 4 chars
        available_data.year = [available_data.year; year];

        % Parse station ID from filename and add to output structure
        station_id = extractBefore( ...
            filename_base, ...
            strlength(filename_base) - 5 + 1);      % grab all but last 5 chars
        available_data.Station_id = [available_data.Station_id; ...
                                     convertCharsToStrings(station_id{1})];
    end
end

% Filter results to a single buoy if a buoy is specified and it exists
if options.buoy_number ~= ""
    indices = ismember(available_data.Station_id, options.buoy_number);
    if any(indices)
        for field_name = fieldnames(available_data)'
            if length(available_data.(field_name{1})) == length(indices)
                available_data.(field_name{1}) = ...
                    available_data.(field_name{1})(indices);
            end
        end
    else
        warning('MATLAB:NDBC_available_data', ...
                'Buoy number "%s" could not be found.', ...
                options.buoy_number);
    end
end
