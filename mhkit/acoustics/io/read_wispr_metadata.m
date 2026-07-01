function metadata = read_wispr_metadata(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads the metadata from the WISPR .dat file and
% returns the metadata in a structure.
%
% Parameters
% ------------
%   filename: string
%       Path to WISPR .dat file.
%
% Returns
% ---------
%   metadata: struct
%       A structure containing metadata such as sampling_rate,
%       adc_vref, gain, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
end

arguments (Output)
    metadata struct
end

% Open WISPR .dat file
fid = fopen(filename, 'rb');
if fid == -1
    error('MHKiT:acoustics:read_wispr_metadata:FileNotFound', 'Could not open file: %s', filename);
end

% Read the 512-byte header
header_raw = fread(fid, 512, 'uint8=>uint8');
fclose(fid);

% Convert raw bytes to standard character array
header_str = char(header_raw');

% Split header string into lines
lines = strsplit(header_str, {'\r\n', '\n', '\r'});

metadata = struct();

for i = 1:length(lines)
    line = strtrim(lines{i});
    
    if isempty(line)
        continue;
    end
    
    % If we hit non-printable binary garbage or null-padding, stop parsing.
    char_vals = double(line);
    if any(char_vals < 32 & char_vals ~= 9) || any(char_vals > 126)
        break;
    end
    
    parts = strsplit(line, '=');
    if length(parts) == 2
        key = strtrim(parts{1});
        value = strtrim(parts{2});
        
        % Strip trailing semicolons
        value = regexprep(value, ';+$', '');
        
        % Normalize key to a valid MATLAB struct field name
        valid_key = matlab.lang.makeValidName(key);
        
        if contains(value, "'")
            % String value, remove single quotes
            value = strrep(value, "'", "");
            metadata.(valid_key) = value;
        else
            % Numeric value
            val_num = str2double(value);
            if ~isnan(val_num)
                metadata.(valid_key) = val_num;
            else
                metadata.(valid_key) = value;
            end
        end
    elseif contains(line, 'WISPR')
        parts_space = strsplit(line, ' ');
        metadata.version = parts_space{end};
    end
end

% Calculate file length if not explicitly defined but other components are present
if ~isfield(metadata, 'file_length_sec') && isfield(metadata, 'file_size') && ...
        isfield(metadata, 'sample_size') && isfield(metadata, 'sampling_rate')
    metadata.file_length_sec = metadata.file_size * 512 / metadata.sample_size / metadata.sampling_rate;
end

end
