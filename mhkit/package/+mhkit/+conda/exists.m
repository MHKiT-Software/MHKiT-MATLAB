function result = exists()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Checks if Conda is installed on the current user's computer
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     result: logical
%         true if Conda is installed and accessible
%         false if Conda is not found
%
% Examples
% ---------
%     is_conda_installed = exists()
%     if is_conda_installed
%         disp('Conda is installed!');
%     else
%         disp('Conda is not installed.');
%     end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % On Windows, try to add conda to PATH first
    if ispc
        userprofile = getenv('USERPROFILE');
        conda_script_paths = {
            fullfile(userprofile, 'miniconda3', 'Scripts')
            fullfile(userprofile, 'miniconda3', 'condabin')
            fullfile(userprofile, 'Miniconda3', 'Scripts') 
            fullfile(userprofile, 'Miniconda3', 'condabin')
            fullfile(userprofile, 'anaconda3', 'Scripts')
            fullfile(userprofile, 'anaconda3', 'condabin')
        };
        
        current_path = getenv('PATH');
        for i = 1:length(conda_script_paths)
            if exist(conda_script_paths{i}, 'dir')
                new_path = [conda_script_paths{i} ';' current_path];
                setenv('PATH', new_path);
                break; % Only add the first found path
            end
        end
    end

    % First, try to check using system path
    try
        % Attempt to run conda info with system command
        [status, out] = mhkit.sys('conda info --base');

        % If the command succeeds (status = 0), Conda is installed
        if status == 0
            result = true;
            return;
        end
    catch
        % If system command fails, continue to next check
    end

    % Check common Conda installation paths
    % List of potential Conda installation directories
    conda_paths = {
        fullfile(getenv('HOME'), 'anaconda3', 'bin', 'conda')
        fullfile(getenv('HOME'), 'miniconda3', 'bin', 'conda')
        fullfile(getenv('USERPROFILE'), 'Anaconda3', 'Scripts', 'conda.exe')
        fullfile(getenv('USERPROFILE'), 'Miniconda3', 'Scripts', 'conda.exe')
        '/opt/anaconda3/bin/conda'
        '/usr/local/anaconda3/bin/conda'
    };

    % Check each potential path
    for i = 1:length(conda_paths)
        if exist(conda_paths{i}, 'file')
            result = true;
            disp(conda_paths{i});
            return;
        end
    end

    % If no paths found, Conda is not installed
    result = false;
end
