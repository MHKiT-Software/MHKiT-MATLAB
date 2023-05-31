function output_name=check_name(input_name)
    output_name = input_name;
    if isempty(regexp(input_name,"^[a-zA-Z]\w*$",'ONCE'))
        fprintf('field name %s is invalid, ',input_name);
        output_name = change_name(input_name);
        fprintf('changed to %s \n',output_name);
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
