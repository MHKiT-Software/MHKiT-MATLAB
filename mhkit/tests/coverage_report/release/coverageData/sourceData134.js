var sourceData134 = {"FileName":"/Users/asimms/Desktop/Programming/mhkit_matlab_simms_dev/MHKiT-MATLAB/mhkit/utils/check_name.m","RawFileContents":["function output_name=check_name(input_name)","%%%%%%%%%%%%%%%%%%%%","%   Check if a given string is a valid","%       field name for MATLAB struct.","%   If not, convert it to a valid field name","%","%   \"Field names can contain ASCII letters (A–Z, a–z), digits (0–9),","%   and underscores, and must begin with a letter.","%   The maximum length of a field name is $namelengthmax$.\"","%   -- MATLAB ref for struct","%","% Parameters","% ------------","%   input_name: string","%       Name to check.","%","% Returns","% ------------","%   output_name: the same as input_name if it is valid,","%       otherwise, it is the modified name.","%   1. All special characters are replaced by '_SC'","%   2. Names not beigining with letters now begins with 'A'.","%   3. Name length exceeding $namelengthmax$ (=63) are truncated","%       to name with a length of $namelengthmax$.","%","%","%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","    output_name = input_name;","    % check if the name is valid:","    % 1. begin with letter; 2. contains only words; 3. <63","    if isempty(regexp(input_name,\"^[a-zA-Z]\\w*$\",'ONCE'))","        fprintf('field name %s is invalid, ',input_name);","        output_name = change_name(input_name);","        fprintf('changed to %s \\n',output_name);","    end","    if strlength(output_name)>namelengthmax","        temp = char(output_name);","        output_name=temp(1:namelengthmax);","    end","    function new_name = change_name(old_name)","        % replace special characters with '_SC'","        new_name = regexprep(old_name,'[^a-zA-Z0-9_]','_SC');","        % make it begin with letters","        if ~isletter(new_name(1))","            new_name = append('A',new_name);","        end","    end","end","",""],"CoverageDisplayDataPerLine":{"Function":[{"LineNumber":1,"Hits":2245,"StartColumnNumbers":0,"EndColumnNumbers":43,"ContinuedLine":false},{"LineNumber":41,"Hits":8,"StartColumnNumbers":4,"EndColumnNumbers":45,"ContinuedLine":false}],"Statement":[{"LineNumber":29,"Hits":2245,"StartColumnNumbers":4,"EndColumnNumbers":29,"ContinuedLine":false},{"LineNumber":32,"Hits":2245,"StartColumnNumbers":4,"EndColumnNumbers":57,"ContinuedLine":false},{"LineNumber":33,"Hits":8,"StartColumnNumbers":8,"EndColumnNumbers":57,"ContinuedLine":false},{"LineNumber":34,"Hits":8,"StartColumnNumbers":8,"EndColumnNumbers":46,"ContinuedLine":false},{"LineNumber":35,"Hits":8,"StartColumnNumbers":8,"EndColumnNumbers":48,"ContinuedLine":false},{"LineNumber":37,"Hits":2245,"StartColumnNumbers":4,"EndColumnNumbers":43,"ContinuedLine":false},{"LineNumber":38,"Hits":3,"StartColumnNumbers":8,"EndColumnNumbers":33,"ContinuedLine":false},{"LineNumber":39,"Hits":3,"StartColumnNumbers":8,"EndColumnNumbers":42,"ContinuedLine":false},{"LineNumber":43,"Hits":8,"StartColumnNumbers":8,"EndColumnNumbers":61,"ContinuedLine":false},{"LineNumber":45,"Hits":8,"StartColumnNumbers":8,"EndColumnNumbers":33,"ContinuedLine":false},{"LineNumber":46,"Hits":0,"StartColumnNumbers":12,"EndColumnNumbers":44,"ContinuedLine":false}]}}