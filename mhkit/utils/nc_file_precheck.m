function nc_file_precheck(filename)
%%%%%%%%%%%%%%%%%%%%
%     Check NetCDF filename before working on it 
%     
% Parameters
% ------------
%   filename: string
%       Filename of NetCDF file to read.
% Returns
% ---------
%   No return.
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check to see if the filename input is a string
    if ~ischar(filename) && ~isstring(filename)
        ME = MException('MATLAB:nc_file_precheck',['filename must be a ' ...
            'character string']);
        throw(ME);
    % check to see if the file exists
    elseif ~isfile(filename)
        ME = MException('MATLAB:nc_file_precheck','file does not exist');
        throw(ME);
    % check MATLAB version & if it is h5 file
    elseif isMATLABReleaseOlderThan("R2021b") || ...
            endsWith(filename, ".h5")
        %ds = read_h5(filename); 
        % return
        ME = MException('MATLAB:nc_file_precheck','please use read_h5');
        throw(ME);
    end

end