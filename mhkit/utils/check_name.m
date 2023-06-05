function output_name=check_name(input_name)
%%%%%%%%%%%%%%%%%%%%
%   Check if a given string is a valid 
%       field name for MATLAB struct. 
%   If not, convert it to a valid field name
%   
%   "Field names can contain ASCII letters (A–Z, a–z), digits (0–9), 
%   and underscores, and must begin with a letter. 
%   The maximum length of a field name is $namelengthmax$." 
%   -- MATLAB ref for struct
%     
% Parameters
% ------------
%   input_name: string
%       Name to check. 
%
% Returns
% ------------
%   output_name: the same as input_name if it is valid, 
%       otherwise, it is the modified name.
%   1. All special characters are replaced by '_SC'
%   2. Names not beigining with letters now begins with 'A'.
%   3. Name length exceeding $namelengthmax$ (=63) are truncated 
%       to name with a length of $namelengthmax$.
% 
%     
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_name = input_name;
    % check if the name is valid:
    % 1. begin with letter; 2. contains only words; 3. <63
    if isempty(regexp(input_name,"^[a-zA-Z]\w*$",'ONCE'))
        fprintf('field name %s is invalid, ',input_name);
        output_name = change_name(input_name);
        fprintf('changed to %s \n',output_name);
    end
    if strlength(output_name)>namelengthmax
        temp = char(output_name);
        output_name=temp(1:namelengthmax);
    end
    function new_name = change_name(old_name)
        % replace special characters with '_SC'
        new_name = regexprep(old_name,'[^a-zA-Z0-9_]','_SC');
        % make it begin with letters
        if ~isletter(new_name(1))
            new_name = append('A',new_name);
        end
    end
end
