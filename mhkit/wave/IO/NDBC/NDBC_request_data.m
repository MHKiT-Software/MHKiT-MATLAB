function ndbc_data=NDBC_request_data(parameter,filenames,options)

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
%     proxy : structure or None (optional)
%         Parameter is now deprecated. Will ignore and print a warning.
%         To request data from behind a firewall, configure in MATLAB Preferences.
%         for example "localhost:8080"
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
    options.proxy string = "";
end

if (options.proxy ~= "")
    warning('To use a proxy server, configure in MATLAB Preferences.');
end

for i = 1:length(filenames)
    % Formulate query
    data_url = "https://www.ndbc.noaa.gov/data/historical/";
    api_query = strcat(parameter, '/', filenames(i));

    % Submit query and get data
    response = webread(data_url + api_query);

    % Write compressed data to temporary file
    compressed_filename = fullfile(tempdir, filenames(i));
    fileID = fopen(compressed_filename, 'w');
    fwrite(fileID, response);
    fclose(fileID);

    % Unzip temporary GNU zip (.gz) file
    uncompressed_filename = gunzip(compressed_filename);

    % Read uncompressed delimited data file
    opts = detectImportOptions(uncompressed_filename{1});
    file_data = readmatrix(uncompressed_filename{1}, opts);

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

    % Parse data and add to output structure
    data_struct.YYYY = file_data(2:end, 1);                     % years
    data_struct.MM = file_data(2:end, 2);                       % months
    data_struct.DD = file_data(2:end, 3);                       % days
    data_struct.hh = file_data(2:end, 4);                       % hours
    data_struct.time = datetime( ...
        year, ...
        data_struct.MM, ...
        data_struct.DD, ...
        data_struct.hh, ...
        0, ...
        0);
    data_struct.frequency = transpose(file_data(1, 5:end));     % column vector
    data_struct.spectrum = transpose(file_data(2:end, 5:end));  % [freq, time]
    ndbc_data.(strcat('ID_', station_id)). ...
        (strcat('year_', num2str(year))) = data_struct;
end

% Remove top structure level if there's only one buoy
field_names = fieldnames(ndbc_data);
if length(field_names) == 1
    ndbc_data = ndbc_data.(field_names{1});
end
