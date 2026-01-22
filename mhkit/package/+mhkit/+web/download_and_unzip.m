function [status, extracted_path] = download_and_unzip(url, target_dir, unzip_folder_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Downloads a zip file from a URL and extracts its contents
%
% Parameters
% ------------
%     url: string
%         URL of the zip file to download
%
%     target_dir: string
%         Directory where the zip file should be downloaded to
%         Will be created if it doesn't exist
%
%     unzip_folder_name: string
%         Name of folder to extract to
%         This folder will be deleted if it exists before extraction
%
% Returns
% ---------
%     status: logical
%         true for success, false for failure
%
%     extracted_path: string
%         Full path to the extracted directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check that all required parameters are provided
    if nargin < 3
        error('All parameters are required: url, target_dir, and unzip_folder_name');
    end

    % Initialize return values
    status = false;
    extracted_path = '';

    try
        % Create target directory if it doesn't exist
        if ~exist(target_dir, 'dir')
            mkdir(target_dir);
            fprintf('Created directory: %s\n', target_dir);
        end

        % Create path for the specified extraction folder
        extraction_dir = fullfile(target_dir, unzip_folder_name);

        % Delete extraction folder if it already exists
        if exist(extraction_dir, 'dir')
            fprintf('Deleting existing folder: %s\n', extraction_dir);
            rmdir(extraction_dir, 's');  % 's' flag for recursive deletion
        end

        % Create the extraction folder
        mkdir(extraction_dir);
        fprintf('Created extraction folder: %s\n', extraction_dir);

        % Generate a temporary filename for the zip file
        zip_file = fullfile(target_dir, 'temp_download.zip');

        % Download the zip file
        fprintf('Downloading %s to %s...\n', url, zip_file);
        websave(zip_file, url);

        % Check if download was successful
        if ~exist(zip_file, 'file')
            fprintf('Failed to download the file.\n');
            return;
        end

        % Unzip the file
        fprintf('Extracting contents to %s...\n', extraction_dir);
        unzip(zip_file, extraction_dir);

        % Set the extracted path to the specified folder
        extracted_path = extraction_dir;

        % Delete the temporary zip file
        delete(zip_file);
        fprintf('Deleted temporary zip file.\n');

        % Set success status
        status = true;
        fprintf('Successfully downloaded and extracted zip file to %s.\n', extracted_path);

    catch ME
        % Handle any errors
        status = false;
        fprintf('Error: %s\n', ME.message);
    end
end
