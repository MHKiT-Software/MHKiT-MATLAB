function ndbc_data=NDBC_request_data(parameter,filenames)

%%%%%%%%%%%%%%%%%%%%
%     Requests data by filenames and returns a structure of structures
%     for each filename passed. If filenames for a single buoy are passed
%     then the yearly structures in the returned structure (ndbc_data) are
%     indexed by year (e.g. ndbc_data.year_2004). If multiple buoy ids are
%     passed then the returned dictionary is indexed by buoy id and year
%     (e.g. ndbc_data['46022']['2014']).
%
%
% Parameters
% ------------
%     parameter : string
%         'swden'  : 'Raw Spectral Wave Current Year Historical Data'
%         'stdmet' : 'Standard Meteorological Current Year Historical Data'
%
%     filenames : array of strings
%         Data filenames on https://www.ndbc.noaa.gov/data/historical/{parameter}/
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
%     ndbc_data: Structure
%         Structure of structures broken down by buoy and years of data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    parameter string {mustBeMember(parameter,{'swden','stdmet'})}
    filenames (1,:) string
end

MAX_RETRIES = 5;                         % number of query retries if error

for i = 1:length(filenames)
    % Formulate query
    data_url = "https://www.ndbc.noaa.gov/data/historical/";
    api_query = strcat(parameter, '/', filenames(i));

    % Submit query and get data
    for j = 0:MAX_RETRIES
        try
            response = webread(data_url + api_query);
            break;
        catch ME
            if i == MAX_RETRIES
                disp(['MATLAB:NDBC_request_data: ', ME.identifier]);
                rethrow(ME)
            else
                pause(1);   % pause(seconds) and retry query
            end
        end
    end

    % Write compressed data to temporary file
    compressed_filename = fullfile(tempdir, filenames(i));
    fileID = fopen(compressed_filename, 'w');
    fwrite(fileID, response);
    fclose(fileID);

    % Unzip temporary GNU zip (.gz) file
    uncompressed_filename = gunzip(compressed_filename);

    % Read uncompressed delimited data file to a table
    % (suppress warning of '.'s being converted to 'x_'s in table header)
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    file_table = readtable(uncompressed_filename{1}, ...
        'ReadVariableNames', true, ...
        'VariableNamesLine', 1, ...
        'ExtraColumnsRule', 'ignore');  % ignore extra columns data w/o headings
    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
    file_header = file_table.Properties.VariableDescriptions;   % includes '.'s

    % Delete compressed and uncompressed files
    delete(compressed_filename);
    delete(uncompressed_filename{1});

    % Parse year and station_id from filename
    [~, filename_no_extension, ~] = fileparts(uncompressed_filename{1});
    if contains(filename_no_extension, '_')
        % Handles filenames like: 42002w2008_old.txt.gz
        filename_base = extractBefore(filename_no_extension, '_');
    else
        filename_base = filename_no_extension;
    end
    year = str2double(extractAfter( ...
        filename_base, ...
        strlength(filename_base) - 4));     % grab last 4 chars
    station_id = extractBefore( ...
        filename_base, ...
        strlength(filename_base) - 5 + 1);  % grab all but last 5 chars

    % Check for bad file data
    if width(file_table) <= 1
        warning('Bad file data for %s, skipping.', filename_base);
        continue;
    end

    % Rename {'YY', '#YY'} columns, handle minutes and add datetime column
    headings = file_table.Properties.VariableNames;
    column_year = find(contains(headings, 'YY'), 1, 'first');
    file_table.Properties.VariableNames{column_year} = 'YYYY';
    column_minutes = find(contains(headings, 'mm'), 1, 'first');
    if isempty(column_minutes)
        minutes = 0;
    else
        minutes = file_table{:, column_minutes};
    end
    time = datetime( ...
        year, ...
        file_table.('MM'), ...
        file_table.('DD'), ...
        file_table.('hh'), ...
        minutes, ...
        0);
    column_first_freq = find(contains(file_header, '.'), 1, 'first');
    file_table = addvars(file_table, time, 'Before', column_first_freq);
    column_first_freq = column_first_freq + 1;

    % Create data structure from the table
    if parameter == "swden"
        % Carry over time columns and add frequency and spectum
        data_struct = table2struct( ...
            file_table(:, 1:(column_first_freq-1)), ...
            'ToScalar', true);
        data_struct.frequency = transpose(cellfun( ...
            @str2num, ...
            file_header((column_first_freq - 1):end))); % note no time column
        data_struct.spectrum = transpose(table2array( ...
            file_table(:, column_first_freq:width(file_table)))); % [freq, time]
    else
        % Keep all columns as-is
        data_struct = table2struct(file_table, 'ToScalar', true);
    end

    % Add data structure to aggregated output structure
    ndbc_data.(strcat('ID_', station_id)). ...
        (strcat('year_', num2str(year))) = data_struct;
end

% Remove top structure level if there's only one buoy
field_names = fieldnames(ndbc_data);
if length(field_names) == 1
    ndbc_data = ndbc_data.(field_names{1});
end

