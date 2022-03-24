function ds=dolfyn_read(filename,options)

%%%%%%%%%%%%%%%%%%%%
%     Read a binary Nortek (e.g., .VEC, .wpr, .ad2cp, etc.) or RDI
%    (.000, .PD0, .ENX, etc.) data file.
%     
% Parameters
% ------------
%     filename: string
%         Filename of instrument file to read.
%     userdata: bool or string (optional)
%         true, false, or string of userdata.json filename (default true)
%         Whether to read the '<base-filename>.userdata.json' file.
%     nens: nan, int, or 2-element array (optional)
%         nan (default: read entire file), int, or 2-element tuple 
%         (start, stop) Number of pings to read from the file.
%
%     call with options -> dolfyn_read(filename, userdata=false, nens=12) 
%
% Returns
% ---------
%     ds: structure 
%         Structure from the binary instrument data
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    arguments
        filename 
        options.userdata = true;
        options.nens = nan;
    end
    
    % check to see if the filename input is a string
    if ~ischar(filename)
        ME = MException('MATLAB:dolfyn_read','filename must be a string');
        throw(ME);
    end
    
    % check to see if the file exists
    if ~isfile(filename)
        ME = MException('MATLAB:dolfyn_read','file does not exist');
        throw(ME);
    end
    
    % check to make sure userdata is bool or string
    if ~isa(options.userdata, 'logical') && ~isa(options.userdata, 'string')
        ME = MException('MATLAB:dolfyn_read','userdata must be a logical or string');
        throw(ME);
    end
    
    % check to make sure nens is numeric or nan
    if ~all(isa(options.nens, 'numeric'))
        ME = MException('MATLAB:dolfyn_read','nens must be numeric or nan');
        throw(ME);
    end
    
    % check to make sure if nens is numeric that its length is equal to 1 or 2   
    if ~isnan(options.nens)
        if length(options.nens) < 1 || length(options.nens) > 2
            ME = MException('MATLAB:dolfyn_read','nens must be a single value or tuple');
            throw(ME);
        end        
    end

    % Loop over binary readers until we find one that works.
    fun_map = struct(  'nortek',{@read_nortek},...                       
                       'signature',{@read_signature},...
                       'rdi',{@read_rdi});

    reader = get_filetype();

    if isnan(reader)
        ME = MException('MATLAB:dolfyn_read',['File %s is not recognized' ...
            [' as a file-type that is readable by DOLfYN. If you think' ...
            ' it should be readable, try using the appropriate read' ...
            ' function (`read_rdi`, `read_nortek`, or `read_signature`) ' ...
            'found in dolfyn.io'], filename]);
        throw(ME);
    else        
        try
            ds = feval(fun_map.(reader),...
                filename,userdata=options.userdata,nens=options.nens);
        catch e
            fprintf("\nError in read: %s\n", e.message);
        end
    end

    function type = get_filetype()
    % Detects whether the file is a Nortek, Signature (Nortek), or RDI
    % file by reading the first few bytes of the file.
    % 
    % Returns
    % =======
    %    NaN - Doesn't match any known pattern
    %    'signature' - for Nortek signature files
    %    'nortek' - for Nortek (Vec, AWAC) files
    %    'RDI' - for RDI files
        
        fid = fopen(filename,'r', 'n', 'UTF-8');      % open disk file    
        bytes = fread(fid, 2, "uint8");
        code = dec2hex(bytes, 2);
        code = convertCharsToStrings(strcat(code(1,:),code(2,:)));
        if strcmpi(code,'7f79') || strcmpi(code,'7f7f')
            type = 'rdi';
        elseif strcmpi(code, 'a50a')
            type = 'signature';
        elseif strcmpi(code, 'a505')
            type = 'nortek';
        else
            type = nan;
        end
    end

end

