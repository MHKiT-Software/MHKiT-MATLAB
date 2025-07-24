function valid_fieldname = sanitize_fieldname(input_str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert a string into a valid MATLAB structure fieldname
%
% Parameters
% ------------
%     input_str: string or char
%         String to be converted into a valid fieldname
%
% Returns
% ---------
%     valid_fieldname: char
%         Valid MATLAB structure fieldname
%
% Notes
% ---------
%     Function performs the following conversions:
%     1. Removes invalid characters
%     2. Ensures it starts with a letter
%     3. Ensures it only contains letters, numbers, and underscores
%     4. Converts camelCase to snake_case
%     5. Handles MATLAB keywords by prefixing with 'x_'
%
% Examples
% ---------
%     sanitize_fieldname('My Field (2)') returns 'my_field_2'
%     sanitize_fieldname('123field') returns 'f123field'
%     sanitize_fieldname('$special@chars#') returns 'special_chars'
%     sanitize_fieldname('myCamelCase') returns 'my_camel_case'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input validation
    if ~ischar(input_str) && ~isstring(input_str)
        error('Input must be a character vector or string');
    end

    % Convert to char array if string
    if isstring(input_str)
        input_str = char(input_str);
    end

    % Convert camelCase to space-separated words
    input_str = regexprep(input_str, '([a-z])([A-Z])', '$1 $2');
    % Convert to lowercase for consistency
    input_str = lower(input_str);

    % Replace spaces and special characters with underscores
    % Keep only alphanumeric characters and underscores
    valid_fieldname = regexprep(input_str, '[^a-z0-9]', '_');

    % Remove consecutive underscores
    valid_fieldname = regexprep(valid_fieldname, '_+', '_');

    % Remove leading and trailing underscores
    valid_fieldname = regexprep(valid_fieldname, '^_+|_+$', '');

    % Ensure the fieldname starts with a letter
    if ~isempty(valid_fieldname)
        if ~isletter(valid_fieldname(1))
            valid_fieldname = ['f' valid_fieldname];
        end
    else
        % If string is empty or contains no valid characters
        valid_fieldname = 'field';
    end

    % Handle the case where the result is a MATLAB keyword
    matlab_keywords = {'break', 'case', 'catch', 'classdef', 'continue', 'else', ...
                'elseif', 'end', 'for', 'function', 'global', 'if', ...
                'otherwise', 'parfor', 'persistent', 'return', 'spmd', ...
                'switch', 'try', 'while'};

    if ismember(valid_fieldname, matlab_keywords)
        valid_fieldname = ['x_' valid_fieldname];
    end
end
