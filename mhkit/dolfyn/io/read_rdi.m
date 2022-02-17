function ds = read_rdi(filename,options)
%%%%%%%%%%%%%%%%%%%%
%     Read a TRDI binary data file.
%     
% Parameters
% ------------
%     filename: string
%         Filename of Nortek file to read.
%     userdata: bool or string (optional)
%         true, false, or string of userdata.json filename (default true)
%         Whether to read the '<base-filename>.userdata.json' file.%     
%     nens: nan, int, or 2-element array (optional)
%         nan (default: read entire file), int, or 2-element tuple 
%         (start, stop) Number of pings to read from the file.
%
%     call with options -> read_rdi(filename, userdata=false, do_checksum=true, nens=12) 
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
        options.do_checksum = false;
        options.nens = nan;
    end
end

