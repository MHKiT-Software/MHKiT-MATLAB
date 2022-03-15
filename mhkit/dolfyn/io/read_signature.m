function ds=read_signature(filename,options)

%%%%%%%%%%%%%%%%%%%%
%     Read a Nortek Signature (.ad2cp) datafile
%     
% Parameters
% ------------
%     filename: string
%         Filename of Nortek file to read.
%     userdata: bool or string (optional)
%         true, false, or string of userdata.json filename (default true)
%         Whether to read the '<base-filename>.userdata.json' file.
%     nens: nan, int, or 2-element array (optional)
%         nan (default: read entire file), int, or 2-element tuple 
%         (start, stop) Number of pings to read from the file.
%
%     call with options -> read_signature(filename, userdata=false, nens=12) 
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
        ME = MException('MATLAB:read_signature',['filename must be a' ...
            ' character string']);
        throw(ME);
    end
    
    % check to see if the file exists
    if ~isfile(filename)
        ME = MException('MATLAB:read_signature','file does not exist');
        throw(ME);
    end
    
    % check to make sure userdata is bool or string
    if ~isa(options.userdata, 'logical') && ~isa(options.userdata, 'string')
        ME = MException('MATLAB:read_signature','userdata must be a logical or string');
        throw(ME);
    end
    
    % check to make sure nens is numeric or nan
    if ~all(isa(options.nens, 'numeric'))
        ME = MException('MATLAB:read_signature','nens must be numeric or nan');
        throw(ME);
    end
    
    % check to make sure if nens is numeric that its length is equal to 1 or 2
    nstart = 0;
    if ~isnan(options.nens)
        if length(options.nens) < 1 || length(options.nens) > 2
            ME = MException('MATLAB:read_signature','nens must be a single value or tuple');
            throw(ME);
        end
        if length(options.nens) == 1            
            npings = options.nens(1);
        else
            nstart = options.nens(1);
            npings = options.nens(2);
        end
    else
        npings = options.nens;
    end

    % read the user data if it was requested using the read_userdata
    % function 
    userdata = read_userdata(filename, options.userdata);

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Initialize Routine
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    % Attempt to read the binary file under the Nortek format    
    fid = fopen(filename,'r', 'n', 'UTF-8');      % open disk file
    % read endian to check if this is a nortek file
    endian1 = fread(fid,2,'uint8','b');
    frewind(fid);
    endian2 = fread(fid,2,'uint8','l');
    if all(endian1 == [165;10])
        endian = 'l';
    elseif all(endian2 == [165;10])
        endian = 'b';
    else
        % this is not a nortek file. Move to next reader
        ME = MException('MATLAB:read_nortek',['could not determine the '...
            'endianness of the file. Are you sure this is a Nortek file?']);
        throw(ME);        
    end
    clearvars endian1 endian2;
    frewind(fid);
    % initialize variables
    filesize = find_filesize();
    fclose(fid);
    index = get_index();
    fid = fopen(filename,'r', 'n', 'UTF-8');      % open disk file
    filehead_config = read_filehead_config_string();
    ens_pos = index.pos(boolarray_firstensemble_ping());
    lastblock_iswhole = calc_lastblock_iswhole();
    config = calc_config();
    burst_readers = init_burst_readers();
    unknown_ID_count = struct();
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        End Initialize
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Read File
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ds = readfile(nstart, npings);
    fclose(fid);    % close the file 

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                       Post Process
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ds = sci_data(ds);
    ds = reorg(ds);
    ds = reduce(ds);

    declin = nan;
    fn = fieldnames(userdata);
    for k=1:numel(fn)
        if( contains(fn{k},'dec') )
            declin = userdata.(fn{k});
        else
            ds.attrs.(fn{k}) = userdata.(fn{k});
        end
    end

    ds = create_dataset(ds);
    ds = set_coords(ds,ds.coord_sys);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(ds, 'orientmat')
        if isfield(ds, 'orientation_down')
            od = ds.orientation_down.data;
        else
            od = nan;
        end
        omat = calc_omat(ds.heading.data,...
                         ds.pitch.data,...
                         ds.roll.data,...
                         od);
        ds.orientmat.data = omat;
        ds.orientmat.dims = { 'time', 'earth', 'inst'};
        ds.orientmat.coords.time = ds.time;
        ds.orientmat.coords.earth = {'E' 'N' 'U'};
        ds.orientmat.coords.inst = {'X' 'Y' 'Z'};        
    end
    
    if ~isnan(declin)
        ds = set_declination(ds, declin);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        End Post Process
    %                        End Function
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Read Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function ds = readfile(ens_start, ens_stop)
        % If the lastblock is not whole, we don't read it.
        % If it is, we do (don't subtract 1)
        nens_total = length(ens_pos) - int32(~lastblock_iswhole);
        if isnan(ens_stop) || ens_stop > nens_total
            ens_stop = nens_total;
        end
        ens_start = int32(ens_start);
        ens_stop = int32(ens_stop);
        nens = ens_stop - ens_start;
        outdat = init_data(ens_start, ens_stop);
        outdat.filehead_config = filehead_config;
        fprintf('Reading file %s\n',filename)
        c = 1;
        c26 = 1;
        if ens_start == 0
            ens_start = 1;
        end
        fseek(fid,ens_pos(ens_start),'bof');
        while true
            try 
                hdr = read_hdr();
            catch
                ds = outdat;
                break
            end            
            id = hdr.id;
            ky = join(["id",string(id)],"_");
            if any([21, 23, 24, 28] == id) % vel, bt, vel_b5, echo
                outdat.(ky) = read_burst(outdat.(ky), ky, c);
            elseif id == 26 % alt_raw (altimeter burst)
                if ~isfield(burst_readers.(ky),'nsamp_index')
                    first_pass = true;
                    burst_readers.(ky).nsamp_index = ...
                        find(strcmp(burst_readers.(ky).names,...
                        'altraw_nsamp'));
                    tmp_idx = burst_readers.(ky).nsamp_index;
                    sz_str = fmt_str(burst_readers.(ky).format(1:tmp_idx-1),...
                        burst_readers.(ky).N(1:tmp_idx-1));
                    [fmt, shift] = py_struct_2_bytes_format(sz_str, true);
                    burst_readers.(ky).nsamp_shift = shift;
                else
                    first_pass = false;
                    tmp_idx = burst_readers.(ky).nsamp_index;
                    shift = burst_readers.(ky).nsamp_shift;
                end
                tmp_idx = tmp_idx + 2;
                fseek(fid,shift,0);
                % Now read the num_samples
                sz = fread(fid,1,"uint32",endian);
                fseek(fid,-shift-4,0);
                if first_pass
                    if isempty(burst_readers.(ky).shape{1,tmp_idx})
                         burst_readers.(ky).shape{1,tmp_idx} = [sz,1];
                    else
                        temp = [sz,1];
                        for qq = 1:numel(burst_readers.(ky).shape{1,tmp_idx})
                            temp(end+1) = ...
                                burst_readers.(ky).shape{1,tmp_idx}(qq);
                        end
                        burst_readers.(ky).shape{1,tmp_idx} = temp;
                    end
                    burst_readers.(ky).N{tmp_idx} = sz;
                    sz_str = fmt_str(burst_readers.(ky).format,...
                        burst_readers.(ky).N);
                    [fmt, nbyte] = py_struct_2_bytes_format(sz_str, true);
                    burst_readers.(ky).nbyte = nbyte;
                    %  Initialize the array
                    outdat.(ky).altraw_samp = zeros([sz,1,...
                        length(outdat.(ky).altraw_samp)]);
                else
                    if sz ~= burst_readers.(ky).N{tmp_idx}
                        ME = MException('MATLAB:read_signature:readfile',...
                            ['The number of samples in this Altimeter Raw' ...
                            'burst is different from prior bursts.']);
                        throwAsCaller(ME)
                    end
                end
                outdat.(ky) = read_burst(outdat.(ky), ky, c);
                outdat.(ky).ensemble(c26) = c;
                c26 = c26 + 1;
            elseif any([22, 27, 29, 30, 31, 35, 36] == id) % avg record, 
                % bt record, DVL, alt record, avg alt_raw record, raw echo,
                % raw echo transmit
                warning(['Unhandled ID: 0x:%X (%d)\n    This ID is not yet' ...
                    ' handled by DOLfYN.\n'], id, id)
                fseek(fid,hdr.sz,0);
            elseif id == 160
                % 0xa0 (i.e., 160) is a 'string data record'
                if ~isfield(outdat, ky)
                    outdat.(ky) = struct();
                end
                temp = read_str(hdr.sz);
                s_id = temp.id;
                str = temp.str;
                clear temp;
                outdat.(ky).(join(["id",string(s_id)],"_")) = str;
            else
                if ~isfield(unknown_ID_count,ky)
                    unknown_ID_count.(ky) = 1;
                    fprintf('Unknown ID: 0x%X', id)
                else
                    unknown_ID_count.(ky) = unknown_ID_count.(ky) + 1;
                end
                fseek(fid,hdr.sz,0);
            end

            c = advance_ens_count(c, ens_start, nens_total);
            if c > nens
                ds = outdat;
                break;
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = read_burst(dat, key, ens)
        out = dat;
        dat_array = read(key);
        fields = fieldnames(out);
        qq = 1;
        for i = 1:numel(burst_readers.(key).N)
            l = length(size(out.(fields{i})));
            otherdims = repmat({':'},1,l-1);
            try
                out.(fields{i})(ens,otherdims{:}) = ...
                    dat_array(qq:qq+(burst_readers.(key).N{i}-1));
            catch
                if l < 4
                    out.(fields{i})(:) = ...
                        dat_array(qq:qq+(burst_readers.(key).N{i}-1));
                else
                    dim = length(out.(fields{i})(ens,:,:,1));
                    for j = 1:4                        
                        out.(fields{i})(ens,:,:,j) = ...
                        dat_array(qq+((j-1)*dim):qq+((j)*dim)-1);
                    end
                end
            end
            qq = qq + burst_readers.(key).N{i};
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = read(key)        
        out = zeros([sum([burst_readers.(key).N{:}]),1]);
        qq = 1;
        for i = 1:length(burst_readers.(key).N)
            if ftell(fid) == filesize
                ME = MException('MATLAB:read_signature:readfile',...
                    'End of File.');
                throwAsCaller(ME)
            end
            if strcmp(burst_readers.(key).format{i}, 'B7x')
                out(qq:qq+(burst_readers.(key).N{i}-1)) = ...
                fread(fid,burst_readers.(key).N{i},'uchar',endian);
                qq = qq + burst_readers.(key).N{i};
                fseek(fid,7,0);
            else
                sz_str = fmt_str(burst_readers.(key).format{i},...
                            burst_readers.(key).N{i});
                [fmt, nbyte] = py_struct_2_bytes_format(sz_str, false);
                out(qq:qq+(burst_readers.(key).N{i}-1)) = ...
                    fread(fid,burst_readers.(key).N{i},fmt,endian);
                qq = qq + burst_readers.(key).N{i};
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = get_index()        
        index_file = append(filename , '.index');
        if ~isfile(index_file)
            create_index_slow(index_file, 2^32)
        end        
        f = fopen(index_file,'r', 'n', 'UTF-8');
        fseek(f, 0, "eof");
        size = ftell(f);
        frewind(f);
        temp = fread(f,10,'char',endian);
        if strcmp(convertCharsToStrings(char(temp)), 'Index Ver:') 
            index_ver = fread(f,1,'uchar',endian);
            fseek(f,1,0);
            out = struct('ens',zeros((size-12)/37,1),...
                'hw_ens',zeros((size-12)/37,1),...
                'pos',zeros((size-12)/37,1),...
                'ID',zeros((size-12)/37,1),...
                'config',zeros((size-12)/37,1),...
                'beams_cy',zeros((size-12)/37,1),...
                'blank',zeros((size-12)/37,1),...
                'year',zeros((size-12)/37,1),...
                'month',zeros((size-12)/37,1),...
                'day',zeros((size-12)/37,1),...
                'hour',zeros((size-12)/37,1),...
                'minute',zeros((size-12)/37,1),...
                'second',zeros((size-12)/37,1),...
                'usec100',zeros((size-12)/37,1),...
                'd_ver',zeros((size-12)/37,1));
            out.index_ver = index_ver;
            for i = 1:(size-12)/37
                out.ens(i) = fread(f,1,'uint64');
                out.hw_ens(i) = fread(f,1,'uint32');
                out.pos(i) = fread(f,1,'uint64');
                out.ID(i) = fread(f,1,'uint16');
                out.config(i) = fread(f,1,'uint16');
                out.beams_cy(i) = fread(f,1,'uint16');
                out.blank(i) = fread(f,1,'uint16');
                out.year(i) = fread(f,1,'uint8');
                out.month(i) = fread(f,1,'uint8');
                out.day(i) = fread(f,1,'uint8');
                out.hour(i) = fread(f,1,'uint8');
                out.minute(i) = fread(f,1,'uint8');
                out.second(i) = fread(f,1,'uint8');
                out.usec100(i) = fread(f,1,'uint16');
                out.d_ver(i) = fread(f,1,'uint8');
            end
        else
            % Pre-versioning the index files
            index_ver = nan;
            frewind(f);
            out = struct('ens',zeros(size/32,1),...
                'pos',zeros(size/32,1),...
                'ID',zeros(size/32,1),...
                'config',zeros(size/32,1),...
                'beams_cy',zeros(size/32,1),...
                'blank',zeros(size/32,1),...
                'year',zeros(size/32,1),...
                'month',zeros(size/32,1),...
                'day',zeros(size/32,1),...
                'hour',zeros(size/32,1),...
                'minute',zeros(size/32,1),...
                'second',zeros(size/32,1),...
                'usec100',zeros(size/32,1));
            out.index_ver = index_ver;
            for i = 1:size/32
                out.ens(i) = fread(f,1,'uint64');
                out.pos(i) = fread(f,1,'uint64');
                out.ID(i) = fread(f,1,'uint16');
                out.config(i) = fread(f,1,'uint16');
                out.beams_cy(i) = fread(f,1,'uint16');
                out.blank(i) = fread(f,1,'uint16');
                out.year(i) = fread(f,1,'uint8');
                out.month(i) = fread(f,1,'uint8');
                out.day(i) = fread(f,1,'uint8');
                out.hour(i) = fread(f,1,'uint8');
                out.minute(i) = fread(f,1,'uint8');
                out.second(i) = fread(f,1,'uint8');
                out.usec100(i) = fread(f,1,'uint16');
            end
        end
        fclose(f);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    function out2 = read_filehead_config_string()
        hdr = read_hdr();
        out = struct;
        temp = read_str(hdr.sz);
        s_id = temp.id;
        string = temp.str;
        clear temp;
        splt_str = splitlines(string);
        for i = 1:numel(splt_str)
            if length(split(splt_str(i),',')) < 2
                continue
            end
            temp = split(splt_str(i),',');            
            ky = temp(1);
            val = temp(2:end);
            if isfield(out,ky)
                tmp = out.(ky);
                out.(ky) = strings(numel(tmp)+numel(val),1);
                out.(ky)(1:numel(tmp)) = tmp;
                out.(ky)(numel(tmp)+1:end) = val;              
            else
                out.(ky) = val;
            end
        end
        out2 = struct();
        fields = fieldnames(out);
        for i = 1:numel(fields)
            if startsWith(fields(i),"GET")
                dat = out.(fields{i});
                d = struct();
                for qq = 1:numel(dat)
                    temp = split(dat(qq),'=');            
                    k = temp(1);
                    val = temp(2:end);                   
                    if ~isnan(str2double(val))
                        val = str2double(val);
                    end
                    d.(k) = val;
                end
                out2.(fields{i}(4:end)) = d;
            else
                out2.(fields{i}) = out.(fields{i});
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    function out = read_hdr()
        bytes1 = fread(fid,4,'uchar',endian);
        bytes2 = fread(fid,2,'ushort',endian);
        bytes3 = fread(fid,1,'uchar',endian);
        out = struct();
        out.sync = bytes1(1);
        out.hsz = bytes1(2);
        out.id = bytes1(3);
        out.fam = bytes1(4);
        out.sz = bytes2(1);
        out.cs = bytes2(2);
        out.hcs = bytes3(1);
        if out.sync ~= 165
            ME = MException('MATLAB:read_signature','Out of sync!');
            throwAsCaller(ME)
        end
        fseek(fid,1,0);
    end
    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Helper Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    function outdat = init_data(ens_start, ens_stop)
        outdat = struct();
        nens = int32(ens_stop - ens_start);
        n26 = 0;
        for i = 1:length(index.ID)
            if index.ID(i) == 26 && index.ens(i) >= ens_start && ...
                    index.ens(i) < ens_stop
                n26 = n26 + 1;
            end
        end
        fields = fieldnames(burst_readers);
        for i = 1:numel(fields)
            ky = fields{i};
            if strcmp(ky,'id_26')
                n = n26;
                ens = zeros([n,1]);
            else
                ens = ens_start:1:ens_stop-1;
                n = nens;
            end
            outdat.(ky) = init_burst_data(n, burst_readers.(ky));
            outdat.(ky).ensemble = ens;
            outdat.(ky).units = burst_readers.(ky).units;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = init_burst_data(npings, burst_st)
        out = struct();
        for i = 1 : numel(burst_st.names)
            nm = burst_st.names{i};
            shape = burst_st.shape{i};
            size_ = [npings];
            if ~isempty(shape)
                for qq = 1:numel(shape)
                    size_(end+1) = shape(qq);
                end
            else
                size_(2) = 1;
            end
            out.(nm) = nan(size_);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function create_index_slow(outfile, N_ens)
        fin = fopen(filename, 'r');
        fout = fopen(outfile, 'w');
        fwrite(fout, 'Index Ver:', 'uchar', 'l');
        fwrite(fout, 1, 'ushort', 'l');
        fseek(fin, 0, "eof");
        size = ftell(fin);
        frewind(fin);
        ens = 0;
        N = 0;
        config = 0;
        last_ens = -1;        
        while N < N_ens
            pos = ftell(fin);
            if pos == size                
                break
            end
            dat1 = fread(fin,4,'uchar',endian);
            dat2 = fread(fin,3,'short',endian);
            
            if any([21, 23, 24, 26, 28] == dat1(3))
                d_ver = fread(fin,1,'uchar',endian);
                d_off = fread(fin,1,'uchar',endian);
                config = fread(fin,1,'ushort',endian);
                fseek(fin,4,0);
                yr = fread(fin,1,'uchar', endian);  
                mo = fread(fin,1,'uchar', endian); 
                dy = fread(fin,1,'uchar', endian); 
                h  = fread(fin,1,'uchar', endian);
                m  = fread(fin,1,'uchar', endian);
                s  = fread(fin,1,'uchar', endian);
                u  = fread(fin,1,'ushort',endian);
                fseek(fin,14,0);
                beams_cy = fread(fin,1,'ushort',endian);
                if dat1(3) ~= 23
                    fseek(fin,40,0);
                else
                    fseek(fin,42,0);
                end
                ens = fread(fin,1, 'uint', endian);
                if dat1(3) == 23 || dat1(3) == 28
                    ens = last_ens;
                end
                if last_ens > 0 && last_ens ~= ens
                    N = N + 1;
                end
                fwrite(fout, N, 'uint64', 'l');
                fwrite(fout, ens, 'uint', 'l');
                fwrite(fout, pos, 'uint64', 'l');
                fwrite(fout,[dat1(3), config, beams_cy, 0], 'ushort', 'l');
                fwrite(fout, [yr, mo + 1, dy, h, m, s], 'uchar', 'l');
                fwrite(fout, u, 'ushort', 'l');
                fwrite(fout, d_ver, 'uchar', 'l');
                if dat1(3) ~= 23
                    fseek(fin,dat2(1) - (36 + 40),0);
                else
                    fseek(fin,dat2(1) - (36 + 42),0); 
                end     
                last_ens = ens;
            else
                fseek(fin,dat2(1),0);
            end
        end
        fclose(fin);
        fclose(fout);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    function size = find_filesize()
        fseek(fid, 0, "eof");
        size = ftell(fid);
        frewind(fid);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out_string = read_str(size)
        out_string = struct();
        out_string.id = fread(fid,1,"int8",endian);
        str = fread(fid,size-1,"char",endian);
        str = convertCharsToStrings(char(str));
        out_string.str = str;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    function out = boolarray_firstensemble_ping()
        if all(index.ens == 0) && all(index.hw_ens == 1)
            out = index.ID(1);
        else
            out = false(size(index.ens));
            out(1) = true;
            out(2:end) = diff(index.ens) ~= 0;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = calc_lastblock_iswhole()
        blocksize = unique(diff(ens_pos));
        counts = histcounts(diff(ens_pos), length(blocksize));
        [m, ind] = max(counts);
        standard_blocksize = blocksize(ind);
        out = (filesize - ens_pos(end)) == standard_blocksize;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function config = calc_config()
        ids = unique(index.ID);
        config = struct();
        for qq = 1:numel(ids)
            if ~any([21, 23, 24, 26, 28] == ids(qq))
                continue
            end
            if ids(qq) == 23
                type = 'bt';
            else
                type = 'burst';
            end
            inds = index.ID == ids(qq);
            config_ = index.config(inds);
            beams_cy_ = index.beams_cy(inds);
            if length(unique(config_)) ~= 1
                ME = MException('MATLAB:read_signature',['config are not'...
                    ' identical for id: 0x%X.', ids(qq)]);
                throwAsCaller(ME)
            end
            if length(unique(beams_cy_)) ~= 1
                ME = MException('MATLAB:read_signature',['config are not'...
                    ' identical for id: 0x%X.', ids(qq)]);
                throwAsCaller(ME)
            end
            key_name = join(["id",string(ids(qq))],"_");
            config.(key_name) = headconfig_int2dict(config_(1), type);
            config.(key_name) = beams_cy_int2dict(beams_cy_(1), ids(qq));
            config.(key_name).config_ = config_(1);
            config.(key_name).beams_cy_ = beams_cy_(1);
            config.(key_name).type = type;
            if isfield(config.(key_name), 'cy')
                config.(key_name) = rmfield(config.(key_name),'cy');
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = headconfig_int2dict(val, mode)
        if strcmpi(mode,'burst')
            out = struct();
            out.press_valid = getbit(val, 0);
            out.temp_valid = getbit(val, 1);
            out.compass_valid = getbit(val, 2);
            out.tilt_valid = getbit(val, 3);
            % bit 4 is unused
            out.vel = getbit(val, 5);
            out.amp = getbit(val, 6);
            out.corr = getbit(val, 7);
            out.alt = getbit(val, 8);
            out.alt_raw = getbit(val, 9);
            out.ast = getbit(val, 10);
            out.echo = getbit(val, 11);
            out.ahrs = getbit(val, 12);
            out.p_gd = getbit(val, 13);
            out.std = getbit(val, 14);
            % bit 15 is unused
        else
            out = struct();
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = getbit(val,n)
        out = bitand(bitshift(val,-n),1);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function strct = beams_cy_int2dict(val, id)
        strct = struct();
        cy_opts = {'ENU', 'XYZ', 'beam', nan};
        if id == 28
            strct.n_cells = val;
        else
            strct.n_cells = bitand(val,2^10 -1);
            strct.cy = cy_opts{bitand(bitshift(val,-10),3) + 1};
            strct.n_beams = bitshift(val,-12);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function burst_readers = init_burst_readers()
        burst_readers = struct();
        fields = fieldnames(config);
        for i = 1:numel(fields)
            rdr_id = fields{i};
            if strcmpi(rdr_id,'id_28')
                burst_readers.(rdr_id) = calc_echo_struct( ...
                    config.(rdr_id).config_, config.(rdr_id).n_cells);
            elseif strcmpi(rdr_id,'id_23')
                burst_readers.(rdr_id) = calc_bt_struct( ...
                    config.(rdr_id).config_, config.(rdr_id).n_beams); 
            else
                burst_readers.(rdr_id) = calc_burst_struct(...
                    config.(rdr_id).config_, config.(rdr_id).n_beams,...
                    config.(rdr_id).n_cells); 
            end 
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function calc_echo_struct(cfg, nc)
        out = struct();
        flags = headconfig_int2dict(cfg,'burst');
        dd = get_burst_hdr_struct();
        dd.blank_dist = struct('format', 'H', 'shape', [], 'sci_func',...
                [0.001, 0], 'units', '', 'N', 1);
        if flags.vel || flags.amp || flags.corr || flags.alt || ...
            flags.ast || flags.alt_raw || flags.p_gd || flags.std
            ME = MException('MATLAB:read_signature',['Echosounder'...
                ' ping contains invalid data?']);
            throwAsCaller(ME)
        end
        if flags.echo
            dd.echo = struct('format', 'H', 'shape', [1,nc], 'sci_func',...
                [0.01, 0], 'units', 'dB', 'N', nc);
        end
        if flags.ahrs
            dd.orientmat = struct('format', 'f', 'shape', [1,3,3],...
                'sci_func', nan, 'units', '', 'N', 9);
            dd.quaternion = struct('format', 'f', 'shape', [1,4],...
                'sci_func', nan, 'units', '', 'N', 4);
            dd.angrt = struct('format', 'f', 'shape', [1,3],...
                'sci_func', [pi / 180, 0], 'units', 'rad/s', 'N', 3);
        end
        out = datadef(dd);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function calc_bt_struct(cfg, nb)
        out = struct();
        flags = headconfig_int2dict(cfg,'bt');
        dd = get_bt_hdr_struct();
        if flags.vel
            dd.vel = struct('format', 'i', 'shape', [1, nb], 'sci_func',...
                nan, 'units', 'm/s', 'N', nb);
        end
        if flags.dist
            dd.dist = struct('format', 'i', 'shape', [1, nb], 'sci_func',...
                [0.001, 0], 'units', 'm', 'N', nb);
        end
        if flags.fom
            dd.fom = struct('format', 'H', 'shape', [1, nb], 'sci_func',...
                nan, 'units', '', 'N', nb);
        end
        if flags.ahrs
            dd.orientmat = struct('format', 'f', 'shape', [1,3,3],...
                'sci_func', nan, 'units', '', 'N', 9);
            dd.quaternion = struct('format', 'f', 'shape', [1,4],...
                'sci_func', nan, 'units', '', 'N', 4);
            dd.angrt = struct('format', 'f', 'shape', [1,3],...
                'sci_func', [pi / 180, 0], 'units', 'rad/s', 'N', 3);
        end
        out = datadef(dd);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = calc_burst_struct(cfg, nb, nc)
        out = struct();
        flags = headconfig_int2dict(cfg,'burst');
        dd = get_burst_hdr_struct();
        if flags.echo
            ME = MException('MATLAB:read_signature',['Echosounder'...
                ' data found in velocity ping?']);
            throwAsCaller(ME)
        end
        if flags.vel
            dd.vel = struct('format', 'h', 'shape', [1,nc,nb], 'sci_func',...
                nan, 'units', 'm/s', 'N', nb*nc);
        end
        if flags.amp
            dd.amp = struct('format', 'B', 'shape', [1,nc,nb], 'sci_func',...
                [0.5, 0], 'units', 'dB', 'N', nb*nc);
        end    
        if flags.corr
            dd.corr = struct('format', 'B', 'shape', [1,nc,nb], 'sci_func',...
                nan, 'units', '%', 'N', nb*nc);
        end
        if flags.alt
            dd.alt_dist = struct('format', 'f', 'shape', [], 'sci_func',...
                [1, 0], 'units', 'm', 'N', 1);
            dd.alt_quality = struct('format', 'H', 'shape', [], 'sci_func',...
                [0.01, 0], 'units', 'dB', 'N', 1);
            dd.alt_status = struct('format', 'H', 'shape', [], 'sci_func',...
                nan, 'units', '', 'N', 1);
        end
        if flags.ast
            dd.ast_dist = struct('format', 'f', 'shape', [], 'sci_func',...
                [1, 0], 'units', 'm', 'N', 1);
            dd.ast_quality = struct('format', 'H', 'shape', [], 'sci_func',...
                [0.01, 0], 'units', 'dB', 'N', 1);
            dd.ast_offset_time = struct('format', 'h', 'shape', [],...
                'sci_func', [0.0001, 0], 'units', 's', 'N', 1);
            dd.ast_pressure = struct('format', 'f', 'shape', [], 'sci_func',...
                nan, 'units', 'dbar', 'N', 1);
            dd.ast_spare = struct('format', 'B7x', 'shape', [], 'sci_func',...
                nan, 'units', '', 'N', 1);
        end
        if flags.alt_raw
            dd.altraw_nsamp = struct('format', 'I', 'shape', [],...
                'sci_func', nan, 'units', '', 'N', 1);
            dd.altraw_dsamp = struct('format', 'H', 'shape', [],...
                'sci_func', [0.0001, 0], 'units', 'm', 'N', 1);
            dd.altraw_samp = struct('format', 'h', 'shape', [],...
                'sci_func', nan, 'units', '', 'N', 1);
        end
        if flags.ahrs
            dd.orientmat = struct('format', 'f', 'shape', [1,3,3],...
                'sci_func', nan, 'units', '', 'N', 9);
            dd.quaternion = struct('format', 'f', 'shape', [1,4],...
                'sci_func', nan, 'units', '', 'N', 4);
            dd.angrt = struct('format', 'f', 'shape', [1,3],...
                'sci_func', [pi / 180, 0], 'units', 'rad/s', 'N', 3);
        end
        if flags.p_gd
            dd.percent_good = struct('format', 'B', 'shape', [1,nc],...
                'sci_func', nan, 'units', '%', 'N', nc);
        end
        if flags.std
            dd.pitch_std = struct('format', 'h', 'shape', [],...
                'sci_func', [0.01, 0], 'units', 'deg', 'N', 1);
            dd.roll_std = struct('format', 'h', 'shape', [],...
                'sci_func', [0.01, 0], 'units', 'deg', 'N', 1);
            dd.heading_std = struct('format', 'h', 'shape', [],...
                'sci_func', [0.01, 0], 'units', 'deg', 'N', 1);
            dd.press_std = struct('format', 'h', 'shape', [],...
                'sci_func', [0.1, 0], 'units', 'dbar', 'N', 1);
            dd.std_spare = struct('format', 'H22x', 'shape', [],...
                'sci_func', nan, 'units', '', 'N', 1);
        end
        out = datadef(dd);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = datadef(dd)
        fields = fieldnames(dd);
        out = struct('names',cell(size(numel(fields))),...
            'format',cell(size(numel(fields))),...
            'shape',cell(size(numel(fields))),...
            'sci_func',cell(size(numel(fields))),...
            'units',cell(size(numel(fields))),...
            'N',cell(size(numel(fields))),...
            'nbyte',0);        
        for i = 1:numel(fields)
            out.names{i} = (fields{i});
            out.format{i} = dd.(fields{i}).format;
            out.shape{i} = dd.(fields{i}).shape;
            out.sci_func{i} = dd.(fields{i}).sci_func;
            out.units{i} = dd.(fields{i}).units;
            out.N{i} = dd.(fields{i}).N;
            if dd.(fields{i}).N == 1
                py_char = dd.(fields{i}).format;
            else
                temp = num2str(dd.(fields{i}).N);
                py_char = strcat(temp,dd.(fields{i}).format);
            end
            [format, bytes] = py_struct_2_bytes_format(...
                py_char, true);
            out.nbyte = out.nbyte + bytes;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function dd = get_burst_hdr_struct()
        dd = struct();
        dd.ver = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.DatOffset = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.config = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.SerialNum = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.year = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.month = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.day = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.hour = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.minute = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.second = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.usec100 = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.c_sound = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.1, 0], 'units', 'm/s', 'N', 1); %sci_func -> [scale, offset]
        dd.temp = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg C', 'N', 1); 
        dd.pressure = struct('format', 'I', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'dbar', 'N', 1); 
        dd.heading = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.pitch = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.roll = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.beam_config = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.cell_size = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'm', 'N', 1);
        dd.blank_dist = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'm', 'N', 1);
        dd.nom_corr = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '%', 'N', 1);
        dd.temp_press = struct('format', 'B', 'shape', [], 'sci_func',...
            [0.2,-20], 'units', 'deg C', 'N', 1);
        dd.batt = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.1, 0], 'units', 'V', 'N', 1);
        dd.mag = struct('format', 'h', 'shape', [1,3], 'sci_func',...
            [0.1, 0], 'units', 'uT', 'N', 3);
        dd.accel = struct('format', 'h', 'shape', [1,3], 'sci_func',...
            [1. / 16384 * 9.81, 0], 'units', 'm/s^2', 'N', 3);
        dd.ambig_vel = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'm/s', 'N', 1);
        dd.data_desc = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.xmit_energy = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', 'dB', 'N', 1);
        dd.vel_scale = struct('format', 'b', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.power_level = struct('format', 'b', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.temp_mag = struct('format', 'h', 'shape', [], 'sci_func',...
            nan, 'units', 'deg C', 'N', 1);
        dd.temp_clock = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg C', 'N', 1);
        dd.error = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.status0 = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', 'binary', 'N', 1);
        dd.status = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', 'binary', 'N', 1);
        dd.ensemble_ = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function dd = get_bt_hdr_struct()
        dd = struct();
        dd.ver = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.DatOffset = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.config = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.SerialNum = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.year = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.month = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.day = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.hour = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.minute = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.second = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.usec100 = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.c_sound = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.1, 0], 'units', 'm/s', 'N', 1); %sci_func -> [scale, offset]
        dd.temp = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg C', 'N', 1); 
        dd.pressure = struct('format', 'I', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'dbar', 'N', 1); 
        dd.heading = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.pitch = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.roll = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg', 'N', 1);
        dd.beam_config = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.cell_size = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'm', 'N', 1);
        dd.blank_dist = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'm', 'N', 1);
        dd.nom_corr = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '%', 'N', 1);
        dd.unused = struct('format', 'B', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.batt = struct('format', 'H', 'shape', [], 'sci_func',...
            [0.1, 0], 'units', 'V', 'N', 1);
        dd.mag = struct('format', 'h', 'shape', [1,3], 'sci_func',...
            nan, 'units', 'gauss', 'N', 3);
        dd.accel = struct('format', 'h', 'shape', [1,3], 'sci_func',...
            [1. / 16384 * 9.81, 0], 'units', 'm/s^2', 'N', 3);
        dd.ambig_vel = struct('format', 'I', 'shape', [], 'sci_func',...
            [0.001, 0], 'units', 'm/s', 'N', 1);
        dd.data_desc = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.xmit_energy = struct('format', 'H', 'shape', [], 'sci_func',...
            nan, 'units', 'dB', 'N', 1);
        dd.vel_scale = struct('format', 'b', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.power_level = struct('format', 'b', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.temp_mag = struct('format', 'h', 'shape', [], 'sci_func',...
            nan, 'units', 'deg C', 'N', 1);
        dd.temp_clock = struct('format', 'h', 'shape', [], 'sci_func',...
            [0.01, 0], 'units', 'deg C', 'N', 1);
        dd.error = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
        dd.status = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', 'binary', 'N', 1);
        dd.ensemble_ = struct('format', 'I', 'shape', [], 'sci_func',...
            nan, 'units', '', 'N', 1);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out_chr = fmt_str(fmt_,N_)
        out_chr = '';
        if iscell(N_)
            for i = 1:numel(fmt_)
                if N_{i} == 1
                    out_chr = char(strcat(out_chr,fmt_{i}));
                else
                    temp = num2str(N_{i});
                    py_char = strcat(temp,fmt_{i});
                    out_chr = char(strcat(out_chr,py_char));
                end
            end
        else
            for i = 1:numel(fmt_)
                if N_ == 1
                    out_chr = char(strcat(out_chr,fmt_));
                else
                    temp = num2str(N_);
                    py_char = strcat(temp,fmt_);
                    out_chr = char(strcat(out_chr,py_char));
                end
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = advance_ens_count(c, ens_start, nens_total)
        try
            %  Checks to makes sure we're not already at the end of the
            %  ens_pos array
            posnow = ens_pos(c + ens_start);
        catch
            % We are at the end of the array, set _posnow
            % We use "+1" here because we want the >= in the while
            % loop to fail for this case so that we go ahead and read
            % the next ping without advancing the ens counter.
            posnow = filesize + 1;
        end
        while ftell(fid) >= posnow
            c = c + 1;
            if c + ens_start + 1 >= nens_total
                break
            end
            try                
                posnow = ens_pos(c + ens_start);
            catch
                posnow = filesize + 1;
            end
        end
        out = c;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = sci_data(dat) 
        out = dat;
        fields = fieldnames(burst_readers);
        for i = 1:numel(fields)
            key1 = fields{i};
            f = burst_readers.(key1).names;
            for j = 1:numel(f)
                key2 = f{j};
                if ~isnan(burst_readers.(key1).sci_func{j})
                    sci_func = burst_readers.(key1).sci_func{j};
                    scale = sci_func(1);
                    offset = sci_func(2);
                    out.(key1).(key2)(:) = (out.(key1).(key2)(:)...
                        + offset) * scale;
                end
            end
            if isfield(out.(key1),'vel') && isfield(out.(key1),'vel_scale')
                out.(key1).vel = out.(key1).vel.* 10.^out.(key1).vel_scale;
            end
        end        
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function outdat = reorg(dat)
        % This function grabs the data from the dictionary of data types
        % (organized by ID), and combines them into a single dictionary.
        outdat = struct('data_vars', struct(), 'coords', struct(),...
            'attrs', struct(), 'units', struct(), 'sys', struct(), ...
            'altraw', struct());
        temp = char(dat.filehead_config.ID(1));
        outdat.attrs.filehead_config = dat.filehead_config;
        outdat.attrs.inst_model = temp(6:end-1);       
        outdat.attrs.inst_make = 'Nortek';
        outdat.attrs.inst_type = 'ADCP';
        outdat.attrs.rotate_vars = {'vel'};
        id_tag = struct('id_21', struct('id', 21, 'tag', ''),...
                        'id_23', struct('id', 23, 'tag', '_bt'),...
                        'id_24', struct('id', 24, 'tag', '_b5'),...
                        'id_26', struct('id', 26, 'tag', '_ast'),...
                        'id_28', struct('id', 28, 'tag', '_echo'));
        id_fields = fieldnames(id_tag);
        for i = 1:numel(id_fields)
            key = id_fields{i};
            id = id_tag.(key).id;
            tag = id_tag.(key).tag;
            if id == 24 || id == 26
                collapse_exclude = [1];
            else
                collapse_exclude = [];
            end
            if ~isfield(dat,key)
                continue;
            end

            for qq = 1:numel(dat.(key).units)
                if ~isfield(outdat.units,burst_readers.(key).names{qq})
                    outdat.units.(burst_readers.(key).names{qq}) = ...
                        dat.(key).units{qq};
                end
            end
            outdat.attrs.(strjoin({'burst_config',tag},'')) = ...
                headconfig_int2dict(collapse(dat.(key).config,...
                'config', collapse_exclude),'burst');
            outdat.coords.(strjoin({'time',tag},'')) = ...
                datenum(datestr(datetime(...
                dat.(key).year(:)+ 1900,...
                dat.(key).month,...
                dat.(key).day,...
                dat.(key).hour,...
                dat.(key).minute,...
                dat.(key).second,...
                dat.(key).usec100)));
            temp = beams_cy_int2dict(collapse(dat.(key).beam_config,...
                'beam_config', collapse_exclude),21);
            outdat.attrs.(strjoin({'n_cells',tag},'')) = temp.n_cells;
            outdat.attrs.(strjoin({'coord_sys_axes',tag},'')) = temp.cy;
            outdat.attrs.(strjoin({'n_beams',tag},'')) = temp.n_beams;
            outdat.attrs.(strjoin({'xmit_energy',tag},'')) = ...
                median(dat.(key).xmit_energy);
            outdat.attrs.(strjoin({'ambig_vel',tag},'')) = ...
                median(dat.(key).ambig_vel);
            
            iter_keys = {'SerialNum', 'cell_size', 'blank_dist', ...
                'nom_corr','data_desc', 'vel_scale', 'power_level'};
            for j = 1:numel(iter_keys)
                ky = iter_keys{j};
                outdat.attrs.(strjoin({ky,tag},'')) = collapse(...
                    dat.(key).(ky), collapse_exclude, ky);
            end

            iter_keys = {'c_sound', 'temp', 'pressure', 'heading',...
                'pitch', 'roll','mag', 'accel', 'batt', 'temp_mag',...
                'temp_clock', 'error','status', 'ensemble'};
            for j = 1:numel(iter_keys)
                ky = iter_keys{j};
                outdat.data_vars.(strjoin({ky,tag},'')) = dat.(key).(ky);
                if strcmp(ky,'ensemble')
                    outdat.data_vars.(strjoin({ky,tag},'')) = ...
                        outdat.data_vars.(strjoin({ky,tag},'')) + 1;
                    outdat.units.(strjoin({ky,tag},'')) = '#';
                end                
            end

            iter_keys = {'vel', 'amp', 'corr', 'prcnt_gd', 'echo',...
                'dist', 'orientmat', 'angrt', 'quaternion', ...
                'ast_pressure', 'alt_dist', 'alt_quality', 'alt_status',...
                'ast_dist', 'ast_quality', 'ast_offset_time',...
                'altraw_nsamp', 'altraw_dsamp', 'altraw_samp',...
                'status0', 'fom', 'temp_press', 'press_std',...
                'pitch_std', 'roll_std', 'heading_std'};
            for j = 1:numel(iter_keys)
                ky = iter_keys{j};                
                if isfield(dat.(key),ky)
                    outdat.data_vars.(strjoin({ky,tag},'')) =...
                        dat.(key).(ky);
                end                
            end 
        end

        % Move 'altimeter raw' data to its own down-sampled structure
        if isfield(dat,'id_26')
            fields = fieldnames(outdat.data_vars);
            for i = 1:numel(fields)
                ky = fields{i};
                if endsWith(ky,"_ast")
                    tmp = split(ky,".");
                    grp = tmp{1};
                    if length(tmp) > 1 && ~isfield(outdat.altraw,grp)
                        outdat.altraw.(grp) = struct();
                    end
                    outdat.altraw.(erase(ky,"_ast")) = ...
                        outdat.data_vars.(ky);
                    outdat.data_vars = rmfield(outdat.data_vars,ky);
                end
            end
        end
        
        if strcmpi(outdat.attrs.coord_sys_axes,'XYZ')
            outdat.attrs.coord_sys = 'inst';
        elseif strcmpi(outdat.attrs.coord_sys_axes,'ENU')
            outdat.attrs.coord_sys = 'earth';
        elseif strcmpi(outdat.attrs.coord_sys_axes,'beam')
            outdat.attrs.coord_sys = 'beam';
        end

        tmp = status2data(outdat.data_vars.status);
        if any([0, 1, 2, 3] == tmp.orient_up(1))
            outdat.attrs.orientation = 'horizontal';
        elseif tmp.orient_up(1) == 4
            outdat.attrs.orientation = 'up';
        elseif tmp.orient_up(1) == 5
            outdat.attrs.orientation = 'down';
        elseif tmp.orient_up(1) == 6
            outdat.attrs.orientation = 'AHRS';
        end

        if tmp.auto_orientation(1) == 0
            outdat.attrs.orient_status = 'fixed';
        elseif tmp.auto_orientation(1) == 1
            outdat.attrs.orient_status = 'auto_UD';
        elseif tmp.auto_orientation(1) == 3
            outdat.attrs.orient_status = 'AHRS-3D';
        end

        keys = {'accel', 'angrt', 'mag'};
        for i = 1:3
            ky = keys{i};
            fields = fieldnames(outdat.data_vars);
            for j = 1:numel(fields)
                dky = fields{j};
                if strcmpi(dky,ky) || startsWith(dky,...
                        strjoin({ky,'_'},''))
                    outdat.attrs.rotate_vars{end+1} = dky;
                end
            end
        end

        if isfield(outdat.data_vars,'vel_bt')
            outdat.attrs.rotate_vars{end+1} = 'vel_bt';
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = reduce(dat)
        % This function takes the output from `reorg`, and further 
        % simplifies the data. Mostly this is combining system, 
        % environmental, and orientation data --- from different data 
        % structures within the same ensemble --- by averaging.       

        % Average these fields
        fields = {'c_sound', 'temp', 'pressure','temp_press',...
            'temp_clock', 'temp_mag','batt'};
        for i = 1:numel(fields)
            ky = fields{i};
            dat.data_vars = reduce_by_average(dat.data_vars,...
                ky, strjoin({ky,'_b5'},''));
        end

        % Angle-averaging is treated separately
        fields = {'heading', 'pitch', 'roll'};
        for i = 1:numel(fields)
            ky = fields{i};
            dat.data_vars = reduce_by_average_angle(dat.data_vars,...
                ky, strjoin({ky,'_b5'},''), true);
        end

        tmp = size(dat.data_vars.vel);
        dat.coords.range = (1:tmp(3)) * dat.attrs.cell_size ...
            + dat.attrs.blank_dist;

        if isfield(dat.data_vars,'vel_b5')
            tmp = size(dat.data_vars.vel_b5);
            dat.coords.range_b5 = (1:tmp(3)) * dat.attrs.cell_size_b5 ...
                + dat.attrs.blank_dist_b5;
        end

        if isfield(dat.data_vars,'echo_echo')
            dat.data_vars.echo = dat.data_vars.echo_echo;
            dat.data_vars = rmfield(dat.data_vars,"echo_echo");
            tmp = size(dat.data_vars.echo);
            dat.coords.range_echo =  (1:tmp(end)) * ...
                dat.attrs.cell_size_echo  + dat.attrs.blank_dist_echo;
        end

        if isfield(dat.data_vars,"orientmat")
            dat.attrs.has_imu = true;
            % Signature AHRS rotation matrix returned in "inst->earth"
            % Change to dolfyn's "earth->inst"
            dat.data_vars.orientmat =  ...
                permute(dat.data_vars.orientmat, [1 2 4 3]);
        else
            dat.attrs.has_imu = false;
        end

        dat.attrs.fs = dat.attrs.filehead_config.BURST.SR;
        dat.attrs.filehead_config.BURST = rmfield(...
            dat.attrs.filehead_config.BURST, "SR");
        tmat = dat.attrs.filehead_config.XFBURST;
        dat.attrs.filehead_config = rmfield(...
            dat.attrs.filehead_config, "XFBURST");
        tm = zeros(tmat.ROWS, tmat.COLS);
        for ii = 1:tmat.ROWS
            for qq = 1:tmat.COLS
                tm(ii,qq) = tmat.(...
                    strjoin({'M',num2str(ii),num2str(qq)},''));
            end
        end
        dat.data_vars.beam2inst_orientmat = tm;

        out = dat;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = collapse(vec, name, exclude)
        % Check that the input vector is uniform, then collapse it to a
        % single value, otherwise raise a warning.
        if is_uniform(vec,[])
            out = vec(1);
        elseif is_uniform(vec, exclude)
            out = something;
        else
            warning(['The variable %s is expected to be uniform, but ' ...
                'it is not.\n',name])
            out = vec(1);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = is_uniform(vec, exclude)        
        if ~isempty(exclude)
            vec(exclude) = [];            
        end
        out = all(vec == vec(1));
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = status2data(val)
        out = struct();
        out.wakeup_state = bitindexer(val, [28,32]);
        out.orient_up = bitindexer(val, [25,28]);
        out.auto_orientation = bitindexer(val, [22,25]);
        out.previous_wakeup_state = bitindexer(val, [18,22]);
        out.last_meas_low_voltage_skip = bitindexer(val, [17]);
        out.active_config = bitindexer(val, [16]);
        out.echo_sounder_index = bitindexer(val, [12,16]);
        out.telemetry_data = bitindexer(val, [11]);
        out.boost_running = bitindexer(val, [10]);
        out.echo_sounder_freq_bin = bitindexer(val, [5,10]);
        % 2, 3, 4 unused
        out.bd_scaling = bitindexer(val, [1]);  % if True: cm scaling of blanking dist
        % 0 unused
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = bitindexer(val, slice)
        if length(slice) == 1
            start = slice(1);
            stop = slice(1)+1;
        else
            start = slice(1);
            stop = slice(2);
        end
        mask = 2 ^ (stop - start) - 1;
        out = bitand(bitshift(val,-start),mask); 
        if mask < 2
            out = logical(out);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = reduce_by_average(data, ky0, ky1)
        out = data;
        if isfield(out, ky1)
            tmp = out.(ky1);
            out = rmfield(out, ky1);
            if isfield(out, ky0)
                out.(ky0) = out.(ky0) + tmp;
                out.(ky0) = out.(ky0)./2;
            else
                out.(ky0) = tmp;
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = reduce_by_average_angle(data, ky0, ky1, degrees)
        out = data;
        if degrees
            rad_fact = pi/180.;
        else 
            rad_fact = 1;
        end
        if isfield(out,ky1)
            if isfield(out, ky0)
                out.(ky0) = angle(exp(data.(ky0).*1i.*rad_fact) ...
                    + exp(data.(ky1).*1i.*rad_fact)) ./ rad_fact;
                out = rmfield(out, ky1);
            else
                out.(ky0) = out.(ky1);
                out = rmfield(out, ky1);
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
end