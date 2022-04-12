function [format, bytes] = py_struct_2_bytes_format(py_format, size_only)
    if isstring(py_format)
        py_format = convertStringsToChars(py_format);
    end
    if length(py_format) == 1
        [format, bytes] = assign_value(py_format);
    else
        bytes = 0;
        last_bytes = nan;
        mltply = nan;
        for i = 1:length(py_format)
            chr = py_format(i);
            if ~isstrprop(chr,'alpha')
                if isnan(mltply)
                    mltply = str2double(chr);
                else
                    temp = num2str(mltply);
                    mltply = str2double(strcat(temp,chr));
                end
                continue
            else
                [format, temp] = assign_value(chr);
                if isnan(mltply)
                    bytes = bytes + temp;
                else
                    bytes = bytes + (temp * mltply);
                end
                mltply = nan;
                if isnan(last_bytes)
                    last_bytes = temp;
                else
                    if ~size_only
                        if last_bytes ~= temp                        
                            ME = MException('MATLAB:py_struct_2_bytes_format', ...
                                ['matlab is not able to mix byte sizes.'...
                                ' IE it is possible to mix padding(x) ' ...
                                'and char(c) because they are both of ' ...
                                'size 1 but it is not possible to read '...
                                'a char and double(d) because double is'...
                                ' of size 8.']);
                            throw(ME);
                        end
                    end
                end
            end
        end
    end

    function [fmt, size] = assign_value(input)
        if strcmp(input,'c')
            fmt = 'char*1';
            size = 1;
        elseif strcmp(input,'x')
            fmt = 'char*1';
            size = 1;
        elseif strcmp(input,'b')
            fmt = 'schar';
            size = 1;
        elseif strcmp(input,'B')
            fmt = 'uchar';
            size = 1;
        elseif strcmp(input,'h')
            fmt = 'short';
            size = 2;
        elseif strcmp(input,'H')
            fmt = 'ushort';
            size = 2;
        elseif strcmp(input,'i')
            fmt = 'int32';
            size = 4;
        elseif strcmp(input,'I')
            fmt = 'uint32';
            size = 4;
        elseif strcmp(input,'l')
            fmt = 'long';
            size = 4;
        elseif strcmp(input,'L')
            fmt = 'ulong';
            size = 4;
        elseif strcmp(input,'q')
            fmt = 'int64';
            size = 8;
        elseif strcmp(input,'Q')
            fmt = 'uint64';
            size = 8;
        elseif strcmp(input,'f')
            fmt = 'float';
            size = 4;
        elseif strcmp(input,'d')
            fmt = 'double';
            size = 8;  
        else
            ME = MException('MATLAB:py_struct_2_bytes_format',['Unsupported'...
                ' python struct charcter: %s',input]);
            throwAsCaller(ME)
        end
    end
end

