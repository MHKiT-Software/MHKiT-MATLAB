function ds = read_rdi(filename,options)
%%%%%%%%%%%%%%%%%%%%
%     Read a TRDI binary data file.
%     
% Parameters
% ------------
%     filename: string
%         Filename of TRDI file to read.
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

    % check to see if the filename input is a string
    if ~ischar(filename)
        ME = MException('MATLAB:read_rdi',['filename must be a ' ...
            'character string']);
        throw(ME);
    end
    
    % check to see if the file exists
    if ~isfile(filename)
        ME = MException('MATLAB:read_rdi','file does not exist');
        throw(ME);
    end
    
    % check to make sure userdata is bool or string
    if ~isa(options.userdata, 'logical') && ~isa(options.userdata, 'string')
        ME = MException('MATLAB:read_rdi','userdata must be a logical or string');
        throw(ME);
    end
    
    % check to make sure nens is numeric or nan
    if ~all(isa(options.nens, 'numeric'))
        ME = MException('MATLAB:read_rdi','nens must be numeric or nan');
        throw(ME);
    end
    
    % check to make sure if nens is numeric that its length is equal to 1 or 2
    nstart = 0;
    if ~isnan(options.nens)
        if length(options.nens) < 1 || length(options.nens) > 2
            ME = MException('MATLAB:read_rdi','nens must be a single value or tuple');
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

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                    year    Initialize Routine
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    debug7f79_ = nan;
    nbyte_ = 0;
    extrabytes = 0;
    search_num_ = 30000;
    fixoffset_ = 0;
    source_ = 0;
    winrivprob_ = false;
    cfac_ = 180 / 2 ^ 31;
    fprintf('Reading file %s\n',filename)
    fid = fopen(filename,'r', 'n', 'UTF-8'); % open disk file

    cfg = struct();
    hdr = struct();

    cfg.name = 'wh-adcp' ;
    cfg.sourceprog = 'instrument';
    cfg.prog_ver = 0;

    pos_ = read_hdr();
    read_cfg();

    fseek(fid,pos_,'bof');
    n_avg = 1;
    ensemble = create_ensemble(n_avg, cfg.n_cells);    
    filesize = find_filesize();
    npings = floor(filesize / (hdr.nbyte + 2 + extrabytes));
    vars_read = {'time'};
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        End Initialize
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Read File
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ds = load_data();
    fclose(fid);    % close the file 

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                       Post Process
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>    
    userdata = read_userdata(filename, options.userdata);
    user_fields = fieldnames(userdata);
    for i = 1:numel(user_fields)
        ds.attrs.(user_fields{i}) = userdata.(user_fields{i});
    end

    if isfield(ds.coords, 'time_gps')
        % GPS data not necessarily sampling at the same rate as ADCP DAQ.
        ds = remove_gps_duplicates(ds);
    end

    % Create xarray like dataset from upper level structure
    ds = create_dataset(ds);
    ds = set_coords(ds,ds.coord_sys);
    
    if ~isfield(ds, 'beam2inst_orientmat')
        ds.beam2inst_orientmat.data = calc_beam_orientmat(...
            ds.attrs.beam_angle,...
            strcmp(ds.attrs.beam_pattern,'convex'), true);
        ds.beam2inst_orientmat.dims = {'x_star', 'x'};
        ds.beam2inst_orientmat.coords = struct('x_star',...
            [1,2,3,4], 'x', [1,2,3,4]);
        ds.coords.x_star = [1,2,3,4];
        ds.coords.x = [1,2,3,4];
    end

    if ~isfield(ds, 'orientmat')
        ds.orientmat.data = calc_orientmat(ds);        
        ds.orientmat.coords.time = ds.time;
        ds.orientmat.coords.inst = {'X' 'Y' 'Z'}; 
        ds.orientmat.coords.earth = {'E' 'N' 'U'};
        ds.orientmat.dims = {'time', 'inst', 'earth'};
    end

    % Check magnetic declination if provided via software and/or userdata
    % If magnetic_var_deg is set, this means that the declination is already
    % included in the heading and in the velocity data.
    declin_set = isfield(ds.attrs, 'declination');
    if declin_set
        declin = ds.attrs.declination;
        ds.attrs = rmfield(ds.attrs, 'declination');
    else
        declin = nan;
    end

    if ds.attrs.magnetic_var_deg ~= 0
        ds.attrs.declination = ds.attrs.magnetic_var_deg;
        ds.attrs.declination_in_orientmat = true;
    end

    if ds.attrs.magnetic_var_deg ~= 0 && declin_set
        warning("magnetic_var_deg is set to %.2f degrees in the " + ...
            "binary file %s, AND declination is set in the " + ...
            "userdata.json file. DOLfYN WILL USE THE VALUE of %.2f" + ...
            " degrees in userdata.json. If you want to use the value " + ...
            "in magnetic_var_deg, delete the value from userdata.json " + ...
            "and re-read the file.", ds.attrs.magnetic_var_deg,...
            filename, declin);
        ds.attrs.declination = declin;
    end

    if ~isnan(declin)
        ds = set_declination(ds, declin);
    end
    
    % VMDAS applies gps correction on velocity in .ENX files only
    end_check = split(filename,'.');
    if strcmpi(end_check{end},'ENX')
        ds.attrs.vel_gps_corrected = true;
    else 
        ds.attrs.vel_gps_corrected = false;
    end

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Read Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function dat = load_data()        
        if isnan(options.nens)
            nens_ = floor(npings / n_avg);
            ens_range = 1:nens_;
        elseif length(options.nens) == 2
            nens_ = floor((npings-nstart)/n_avg);
            ens_range = nstart:npings;
            fseek(fid,(hdr.nbytes + 2 + extrabytes)*nstart, 0);
        else
            nens_ = options.nens;
            ens_range = 1:nens_;
        end
        dat = init_data(nens_);
        dat.coords.range = (0:cfg.n_cells-1)...
            *cfg.cell_size + cfg.bin1_dist_m;
        fields = fieldnames(cfg);
        for i = 1:numel(fields)
            nm = fields{i};
            dat.attrs.(nm) = cfg.(nm);
        end

        % Changing each time to epoch time is very slow to do each
        % iteration so we save time seperately and change it at the end.
        time = NaT(size(dat.coords.time));
        ms = zeros(size(dat.coords.time));
       
        for i = 1:nens_
            try
                read_buffer();
            catch
                dat = remove_end(dat, i);
                break;
            end
            temp = ensemble.vel.data(:) == -32.768;
            ensemble.vel.data(temp) = nan;
            % Fix the 'real-time-clock' century
            if ensemble.rtc.data(1) < 100
                ensemble.rtc.data(1) = ensemble.rtc.data(1) + 2000;
            end
            for qq = 1:numel(vars_read)
                nm = vars_read{qq};
                group = ensemble.(nm).group;
                l = length(size((dat.(group).(nm))));
                otherdims = repmat({':'},1,l-1);
                if strcmpi(nm,'time')
                    continue;
                else                    
                    dat.(group).(nm)(i,otherdims{:}) = ...
                        ensemble.(nm).data;
                end
            end
            try 
                clock = ensemble.rtc.data(:);
                clock(7) = clock(7) * 0.01;    
                time(i) = datetime(clock(1)...
                     ,clock(2),clock(3),clock(4),clock(5),clock(6));
                ms(i) = clock(7);
            catch
                time(i) = NaT;
                warning("Invalid time stamp in ping %d.",...
                    ensemble.number);
            end
        end
        % Now convert time to epochtime
        ind = ~isnat(time);
        dat.coords.time(ind) = double(convertTo(time(ind)...
            ,'epochtime','Epoch','1970-01-01'));
        dat.coords.time(ind) = dat.coords.time(ind) + ms(ind);        
        %
        dat = finalize(dat);
        if isfield(dat.data_vars, 'vel_bt')
            dat.attrs.rotate_vars{end+1} = 'vel_bt';
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_buffer()
        ensemble.k = 0;
        while ensemble.k < n_avg
            search_buffer();
            startpos = ftell(fid) - 2;
            read_hdrseg();
            byte_offset = nbyte_ + 2;
            for n = 1:length(hdr.dat_offsets)
                id = fread(fid,1,"uint16");
                winrivprob_ = false;
                retval = read_dat(id);
                if ~retval 
                    break
                end
                byte_offset = byte_offset + nbyte_;
                if n < length(hdr.dat_offsets)
                    oset = hdr.dat_offsets(n+1) - byte_offset;
                    if oset ~= 0
                        fseek(fid,oset,0);
                    end
                    byte_offset = hdr.dat_offsets(n+1);
                else
                    if hdr.nbyte - 2 ~= byte_offset
                        if ~winrivprob_
                            fseek(fid,hdr.nbyte - 2 - byte_offset, 0);
                        end
                    end
                    byte_offset = hdr.nbyte - 2;
                end
            end
            readbytes = ftell(fid) - startpos;
            offset = hdr.nbyte + 2 - byte_offset;
            check_offset(offset);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function retval = read_dat(id)
        retval = true;
        switch id
            case 0
                read_fixed();
            case 128
                read_var();
            case 256
                read_vel();
            case 512
                read_corr();
            case 768
                read_amp();
            case 1024
                read_prcnt_gd();
            case 1280
                read_status();
            case 1536
                read_bottom();
            case 8192
                read_vmdas();
            case 8226
                retval = read_winriver2();
            case 8448
                read_winriver(38);
            case 8449
                read_winriver(97);
            case 8450
                read_winriver(45);
            case 8451
                read_winriver(60);
            case 8452
                read_winriver(38);
            % Loading of these data is currently not implemented:
            case 1793
                skip_Ncol(4);
            case 1794
                skip_Ncol(4);
            case 1795
                skip_Ncol(4);
            case 2560
                skip_Ncol(nan);
            % 0301 Beam 5 Number of good pings
            case 769
                skip_Ncol(nan);
            % 0302 Beam 5 Sum of squared velocities
            case 770
                skip_Ncol(nan);
            % 0303 Beam 5 Sum of velocities
            case 771
                skip_Ncol(nan);
            % 020C Ambient sound profile
            case 524
                skip_Nbyte(4);
            case 12288
                skip_Nbyte(32);
            otherwise
                read_nocode(id);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = read_hdr()
        cfgid = fread(fid,2,'uint8');
        while (cfgid(1) ~= 127 || cfgid(2) ~= 127) || ~checkheader()
            nextbyte  = fread(fid,1,'uint8');
            pos = ftell(fid);
            cfgid(2) = cfgid(1);
            cfgid(1) = nextbyte;
            if rem(pos,1000) == 0
                fprintf(['Still looking for valid cfgid at file' ...
                    ' position %d ...'],pos)
            end
        end
        out = ftell(fid) - 2;
        read_hdrseg();
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_hdrseg()
        hdr.nbyte = fread(fid,1,'int16');
        fseek(fid,1,0);
        ndat = fread(fid,1,'int8');
        hdr.dat_offsets = fread(fid,ndat,'int16');
        nbyte_ = 4 + ndat * 2;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_cfg()
        cfgid = fread(fid,1,'uint16');
        read_cfgseg();
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_cfgseg()
        cfgstart = ftell(fid);
        tmp = fread(fid,5,'uint8');
        cfg.prog_ver = tmp(1) + tmp(2) / 100.;
        switch tmp(1)
            case {4,5}
                cfg.name = 'bb-adcp';
            case {8,9,16}
                cfg.name = 'wh-adcp';
            case {14,23}
                cfg.name = 'os-adcp';
            otherwise
                cfg.name = 'unrecognized firmware version';
        end
        config = tmp(3:4);
        beam_angle_options = [15, 20, 30];
        cfg.beam_angle = beam_angle_options(bitand(config(2),3)+1);
        freq_options = [75, 150, 300, 600, 1200, 2400, 38];
        cfg.freq = freq_options(bitand(config(1),7)+1);
        beam_pattern_options = {'concave','convex'};
        cfg.beam_pattern = beam_pattern_options{...
            int32((bitand(config(1),8)==8)+1)};
        orientation_options = {'down', 'up'};
        cfg.orientation = orientation_options{...
            int32((bitand(config(1),128)==128)+1)};
        fseek(fid,1,0);
        temp = fread(fid,2,'uint8');
        cfg.n_beams = temp(1);
        cfg.n_cells = temp(2);
        temp = fread(fid,3,'uint16');
        cfg.pings_per_ensemble = temp(1);
        cfg.cell_size = temp(2) * .01;
        cfg.blank = temp(3) * .01;
        temp = fread(fid,4,'uint8');
        cfg.prof_mode = temp(1);
        cfg.corr_threshold = temp(2);
        cfg.prof_codereps = temp(3);
        cfg.min_pgood = temp(4);
        cfg.evel_threshold = fread(fid,1,'uint16');
        cfg.sec_between_ping_groups = ...
            sum(fread(fid,3,'uint8').*[60.; 1.; .01]);
        coord_sys = fread(fid,1,'uint8');
        coord_sys_options = {'beam', 'inst', 'ship', 'earth'};
        cfg.coord_sys = coord_sys_options{...
            bitand(bitshift(coord_sys,-3),3)+1};
        cfg.use_pitchroll = bitand(coord_sys, 4) == 4;
        cfg.use_3beam = bitand(coord_sys, 2) == 2; 
        cfg.bin_mapping = bitand(coord_sys, 1) == 1;
        temp = fread(fid,2,'int16');
        cfg.xducer_misalign_deg = temp(1) * .01;
        cfg.magnetic_var_deg = temp(2) * .01;
        temp = fread(fid,2,'uint8');
        cfg.sensors_src = dec2bin(temp(1),8);
        cfg.sensors_avail = dec2bin(temp(2),8);
        temp = fread(fid,2,'uint16');
        cfg.bin1_dist_m = temp(1) * .01;
        cfg.xmit_pulse  = temp(2) * .01;
        temp = fread(fid,3,'uint8');
        cfg.water_ref_cells = temp(1:2);
        cfg.fls_target_threshold = temp(3);
        fseek(fid,1,0);
        cfg.xmit_lag_m = fread(fid,1,'uint16') * .01;
        nbyte_ = 40;
        cfg.configsize = ftell(fid) - cfgstart;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_fixed()
        if isfield(cfg, 'configsize')
            fseek(fid,cfg.configsize, 0);
            nbyte_ = cfg.configsize;
        else
            read_cfgseg();
        end
        nbyte_ = nbyte_ + 2;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_var()
        ensemble.k = ensemble.k + 1;
        vars_to_add = {'number',...
                       'rtc',...
                       'number',...
                       'builtin_test_fail',...
                       'c_sound',...
                       'depth',...
                       'heading',...
                       'pitch',...
                       'roll',...
                       'salinity',...
                       'temp',...
                       'min_preping_wait',...
                       'heading_std',...
                       'pitch_std',...
                       'roll_std',...
                       'adc'};
        for qq = 1:numel(vars_to_add)
            if ~any(strcmp(vars_read, vars_to_add{qq}))
                vars_read{end+1} = vars_to_add{qq};
            end
        end
        ensemble.number.data(ensemble.k) = fread(fid,1,'uint16');
        ensemble.rtc.data(ensemble.k,1,:) = fread(fid,7,'uint8');
        ensemble.number.data(ensemble.k) = ensemble.number.data(ensemble.k)...
            + 65535 * fread(fid,1,'uint8');
        temp = fread(fid,4,'uint16');
        ensemble.builtin_test_fail.data(ensemble.k) = temp(1);
        ensemble.c_sound.data(ensemble.k)   = temp(2);
        ensemble.depth.data(ensemble.k)     = temp(3) * .1;
        ensemble.heading.data(ensemble.k)   = temp(4) * .01;
        temp = fread(fid,4,'int16');
        ensemble.pitch.data(ensemble.k)     = temp(1) * .01;
        ensemble.roll.data(ensemble.k)      = temp(2) * .01;
        ensemble.salinity.data(ensemble.k)  = temp(3);
        ensemble.temp.data(ensemble.k)      = temp(4) * .01;        
        ensemble.min_preping_wait.data(ensemble.k) = ...
            sum(fread(fid,3,'uint8').*[60.; 1.; .01]);
        temp = fread(fid,11,'uint8');
        ensemble.heading_std.data(ensemble.k)   = temp(1);
        ensemble.pitch_std.data(ensemble.k)     = temp(2) * 0.1;
        ensemble.roll_std.data(ensemble.k)      = temp(3) * 0.1;
        ensemble.adc.data(ensemble.k,1,:)       = temp(4:11);
        nbyte_ = 2 + 40;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_vel()
        if ~any(strcmp(vars_read, 'vel'))
            vars_read{end+1} = 'vel';
        end        
        temp = fread(fid,4*cfg.n_cells,'int16') * .001;
        ensemble.vel.data(ensemble.k,1,:,:) = ...
            reshape(temp,[ensemble.k,1,4,cfg.n_cells]);
        nbyte_ = 2 + 4 * cfg.n_cells * 2;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_corr()
        if ~any(strcmp(vars_read, 'corr'))
            vars_read{end+1} = 'corr';
        end 
        temp = fread(fid, 4*cfg.n_cells, 'uint8');
        ensemble.corr.data(ensemble.k,1,:,:) = ...
            reshape(temp,[ensemble.k,1,4,cfg.n_cells]);
        nbyte_ = 2 + 4 * cfg.n_cells;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_amp()
        if ~any(strcmp(vars_read, 'amp'))
            vars_read{end+1} = 'amp';
        end 
        temp = fread(fid, 4*cfg.n_cells, 'uint8');
        ensemble.amp.data(ensemble.k,1,:,:) = ...
            reshape(temp,[ensemble.k,1,4,cfg.n_cells]);
        nbyte_ = 2 + 4 * cfg.n_cells;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_prcnt_gd()
        if ~any(strcmp(vars_read, 'prcnt_gd'))
            vars_read{end+1} = 'prcnt_gd';
        end 
        temp = fread(fid, 4*cfg.n_cells, 'uint8');
        ensemble.prcnt_gd.data(ensemble.k,1,:,:) = ...
            reshape(temp,[ensemble.k,1,4,cfg.n_cells]);
        nbyte_ = 2 + 4 * cfg.n_cells;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_status()
        if ~any(strcmp(vars_read, 'status'))
            vars_read{end+1} = 'status';
        end 
        temp = fread(fid, 4*cfg.n_cells, 'uint8');
        ensemble.status.data(ensemble.k,1,:,:) = ...
            reshape(temp,[ensemble.k,1,4,cfg.n_cells]);
        nbyte_ = 2 + 4 * cfg.n_cells;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_bottom()
        vars_to_add = {'dist_bt', 'vel_bt', 'corr_bt', 'amp_bt',...
            'prcnt_gd_bt'};
        for qq = 1:numel(vars_to_add)
            if ~any(strcmp(vars_read, vars_to_add{qq}))
                vars_read{end+1} = vars_to_add{qq};
            end            
        end
        if source_ == 2
            warning('lat/lon bottom reading passed')
        else
            fseek(fid,14,0);
        end
        ensemble.dist_bt.data(ensemble.k,1,:) = fread(fid,4,'uint16')*0.01;
        ensemble.vel_bt.data(ensemble.k,1,:) = fread(fid,4,'int16')*0.001;
        ensemble.corr_bt.data(ensemble.k,1,:) = fread(fid,4,'uint8');
        ensemble.amp_bt.data(ensemble.k,1,:) = fread(fid,4,'uint8');
        ensemble.prcnt_gd_bt.data(ensemble.k,1,:) = fread(fid,4,'uint8');
        if source_ == 2
            warning('lat/lon bottom reading passed')
        else
            fseek(fid,71 - 45,0);
            nbyte_ = 2 + 68;
        end
        if cfg.prog_ver >= 5.3
            fseek(fid, 78-71, 0);
            ensemble.dist_bt.data(ensemble.k,1,:) = ...
                ensemble.dist_bt.data(ensemble.k,1,:) + ...
                reshape(fread(fid,4,'uint8') * 655.36,...
                size(ensemble.dist_bt.data));
            nbyte_ = nbyte_ + 11;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_vmdas()
        % Read something from VMDAS 
        % The raw files produced by VMDAS contain a binary navigation data
        % block.
        cfg.sourceprog = 'VMDAS';
        source_ = 1;
        vars_to_add = {'time_gps',...
                       'latitude_gps',...
                       'longitude_gps',...
                       'etime_gps',...
                       'elatitude_gps',...
                       'elongitude_gps',...
                       'flags',...
                       'ntime'};
        for qq = 1:numel(vars_to_add)
            if ~any(strcmp(vars_read, vars_to_add{qq}))
                vars_read{end+1} = vars_to_add{qq};
            end 
        end
        utim = fread(fid,4,'uint8');        
        milliseconds = int32(fread(fid,1,'uint32') / 10);
        date = datetime(utim(3) + utim(4) * 256, utim(2), utim(1), ...
            0, 0, milliseconds/1000);
        fseek(fid, 4, 0);
        ensemble.time_gps.data(ensemble.k) = ...
            datenum(datestr(date),'dd-mmm-yyyy HH:MM:SS');
        ensemble.latitude_gps.data(ensemble.k) = fread(fid,1,'int32') ... 
            * cfac_;
        ensemble.longitude_gps.data(ensemble.k) = fread(fid,1,'int32') ... 
            * cfac_;
        milliseconds = int32(fread(fid,1,'uint32') * 10);
        date = datetime(utim(3) + utim(4) * 256, utim(2), utim(1), ...
            0, 0, milliseconds/1000);
        ensemble.etime_gps.data(ensemble.k) = ...
            datenum(datestr(date),'dd-mmm-yyyy HH:MM:SS'); 
        ensemble.elatitude_gps.data(ensemble.k) = fread(fid,1,'int32') ... 
            * cfac_;        
        ensemble.elongitude_gps.data(ensemble.k) = fread(fid,1,'int32') ... 
            * cfac_;        
        fseek(fid, 12, 0);
        ensemble.flags.data(ensemble.k) = fread(fid,1,'uint16');
        fseek(fid, 6, 0);
        utim = fread(fid,4,'uint8');        
        milliseconds = int32(fread(fid,1,'uint32') / 10);
        date = datetime(utim(1) + utim(2) * 256, utim(4), utim(3),...
            0, 0, milliseconds/1000);
        ensemble.ntime.data(ensemble.k) = ...
            datenum(datestr(date),'dd-mmm-yyyy HH:MM:SS'); 
        fseek(fid, 16, 0);
        nbyte_ = 2 + 76;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = read_winriver2()
        out = true;
        startpos = ftell(fid);
        winrivprob_ = true;
        cfg.sourceprog = 'WINRIVER';
        source_ = 3;
        spid = fread(fid, 1, "uint16");
        if spid == 104
            sz = fread(fid, 1, "uint16");
            dtime = fread(fid, 1, "float64");
            start_string = char(fread(fid,6,"char"));
            dummy = char(fread(fid,1,"char"));
            if ~strcmpi('$GPGGA', convertCharsToStrings(start_string))
                out = false;
                return
            end
            gga_time = char(fread(fid,9,"char"));
            time = [str2double(gga_time(1:2)),...
                str2double(gga_time(3:4)),...
                str2double(gga_time(5:6)) + ...
                (str2double(gga_time(8:end))*100)*0.001];
            if ensemble.rtc.data(1) < 100
                ensemble.rtc.data(1) = ensemble.rtc.data(1) + 2000;
            end
            tmp_time = ...
                datetime(cat(2, squeeze(ensemble.rtc.data(1:3))', time)); 
            ensemble.time_gps.data(ensemble.k) = ...
                datenum(datestr(tmp_time),'dd-mmm-yyyy HH:MM:SS');
            fseek(fid, 1, 0);
            ensemble.latitude_gps.data(ensemble.k) = ...
                fread(fid, 1, "float64");
            tcNS = char(fread(fid,1,"char"));
            if strcmp(convertCharsToStrings(tcNS), 'S')
                ensemble.latitude_gps.data(ensemble.k) = ...
                    ensemble.latitude_gps.data(ensemble.k) * -1;
            elseif ~strcmp(convertCharsToStrings(tcNS), 'N')
                out = false;
                return
            end

            ensemble.longitude_gps.data(ensemble.k) = ...
                fread(fid, 1, "float64");
            tcEW = char(fread(fid,1,"char"));
            if strcmp(convertCharsToStrings(tcEW), 'W')
                ensemble.longitude_gps.data(ensemble.k) = ...
                    ensemble.longitude_gps.data(ensemble.k) * -1;
            elseif ~strcmp(convertCharsToStrings(tcEW), 'E')
                out = false;
                return
            end
            [~, ~] = fread(fid, 2, "uint8");
            tmp = fread(fid, 2, "float32");
            ensemble.hdop = tmp(1);
            ensemble.altitude = tmp(2);
            if ~strcmp(char(fread(fid,1,"char")), 'M')
                out = false;
                return
            end
            ggeoid_sep = fread(fid, 1, "float32");
            if ~strcmp(char(fread(fid,1,"char")), 'M')
                out = false;
                return
            end
            gage = fread(fid, 1, "float32");
            gstation_id = fread(fid, 1, "uint16");
            % 4 unknown bytes (2 reserved+2 checksum?)
            % 78 bytes for GPGGA string (including \r\n)
            % 2 reserved + 2 checksum
            vars_to_add = {'longitude_gps', 'latitude_gps', 'time_gps'};
            for qq = 1:numel(vars_to_add)
                if ~any(strcmp(vars_read, vars_to_add{qq}))
                    vars_read{end+1} = vars_to_add{qq};
                end 
            end
            nbyte_ = ftell(fid) - startpos + 2;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_winriver(n)
        winrivprob_ = true;
        cfg.sourceprog = 'WINRIVER';
        if source_ ~= 2 && source_ ~= 3
            source_ = 2;
        end
        startpos = ftell(fid);
        sz = fread(fid, 1, 'uint16');
        tmp = char(fread(fid, sz, "char"));
        nbyte_ = ftell(fid) - startpos + 2;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function read_nocode(id)
        fprintf('Unrecognized ID code: %0.4X.\n', id)
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Helper Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function valid = checkheader()
        valid = false;
        numbytes = fread(fid,1,'int16');
        if numbytes > 0
            fseek(fid, numbytes - 2, 0);
            cfgid = fread(fid,2,'uint8');
            if length(cfgid) == 2
                fseek(fid,-numbytes-2,0);
                if cfgid(1) == 127 && (cfgid(2) == 127 || cfgid(2) == 121)
                    if cfgid(2) == 121 && isnan(debug7f79_)
                        debug7f79_ = true;                        
                    end
                    valid = true;
                end
            end
        else
            fseek(fid, -2, 0);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = create_ensemble(navg, ncells)
        out = struct();

        out.number.data = zeros(navg);
        out.number.dtype = "uint32";
        out.number.group = "data_vars";
        out.number.units = "";

        out.rtc.data = zeros(navg, 1, 7);
        out.rtc.dtype = "uint16";
        out.rtc.group = "sys";
        out.rtc.units = "";

        out.builtin_test_fail.data = false(navg);
        out.builtin_test_fail.dtype = "bool";
        out.builtin_test_fail.group = "data_vars";
        out.builtin_test_fail.units = "";

        out.c_sound.data = zeros(navg);
        out.c_sound.dtype = "float32";
        out.c_sound.group = "data_vars";
        out.c_sound.units = "m/s";

        out.depth.data = zeros(navg);
        out.depth.dtype = "float32";
        out.depth.group = "data_vars";
        out.depth.units = "m";

        out.pitch.data = zeros(navg);
        out.pitch.dtype = "float32";
        out.pitch.group = "data_vars";
        out.pitch.units = "deg";

        out.roll.data = zeros(navg);
        out.roll.dtype = "float32";
        out.roll.group = "data_vars";
        out.roll.units = "deg";

        out.heading.data = zeros(navg);
        out.heading.dtype = "float32";
        out.heading.group = "data_vars";
        out.heading.units = "deg";

        out.temp.data = zeros(navg);
        out.temp.dtype = "float32";
        out.temp.group = "data_vars";
        out.temp.units = "C";

        out.salinity.data = zeros(navg);
        out.salinity.dtype = "float32";
        out.salinity.group = "data_vars";
        out.salinity.units = "psu";

        out.min_preping_wait.data = zeros(navg);
        out.min_preping_wait.dtype = "float32";
        out.min_preping_wait.group = "data_vars";
        out.min_preping_wait.units = "s";

        out.heading_std.data = zeros(navg);
        out.heading_std.dtype = "float32";
        out.heading_std.group = "data_vars";
        out.heading_std.units = "deg";

        out.pitch_std.data = zeros(navg);
        out.pitch_std.dtype = "float32";
        out.pitch_std.group = "data_vars";
        out.pitch_std.units = "deg";

        out.roll_std.data = zeros(navg);
        out.roll_std.dtype = "float32";
        out.roll_std.group = "data_vars";
        out.roll_std.units = "deg";

        out.adc.data = zeros(navg, 1, 8);
        out.adc.dtype = "uint8";
        out.adc.group = "sys";
        out.adc.units = "";

        out.error_status_wd.data = zeros(navg);
        out.error_status_wd.dtype = "float32";
        out.error_status_wd.group = "attrs";
        out.error_status_wd.units = "";

        out.pressure.data = zeros(navg);
        out.pressure.dtype = "float32";
        out.pressure.group = "data_vars";
        out.pressure.units = "dbar";

        out.pressure_std.data = zeros(navg);
        out.pressure_std.dtype = "float32";
        out.pressure_std.group = "data_vars";
        out.pressure_std.units = "dbar";

        out.vel.data = zeros(navg, 1, 4, ncells);
        out.vel.dtype = "float32";
        out.vel.group = "data_vars";
        out.vel.units = "m/s";

        out.amp.data = zeros(navg, 1, 4, ncells);
        out.amp.dtype = "uint8";
        out.amp.group = "data_vars";
        out.amp.units = "counts";

        out.corr.data = zeros(navg, 1, 4, ncells);
        out.corr.dtype = "uint8";
        out.corr.group = "data_vars";
        out.corr.units = "counts";

        out.prcnt_gd.data = zeros(navg, 1, 4, ncells);
        out.prcnt_gd.dtype = "uint8";
        out.prcnt_gd.group = "data_vars";
        out.prcnt_gd.units = "%";

        out.status.data = zeros(navg, 1, 4, ncells);
        out.status.dtype = "float32";
        out.status.group = "data_vars";
        out.status.units = "";

        out.dist_bt.data = zeros(navg, 1, 4);
        out.dist_bt.dtype = "float32";
        out.dist_bt.group = "data_vars";
        out.dist_bt.units = "m";

        out.vel_bt.data = zeros(navg, 1, 4);
        out.vel_bt.dtype = "float32";
        out.vel_bt.group = "data_vars";
        out.vel_bt.units = "m/s";

        out.corr_bt.data = zeros(navg, 1, 4);
        out.corr_bt.dtype = "uint8";
        out.corr_bt.group = "data_vars";
        out.corr_bt.units = "counts";

        out.amp_bt.data = zeros(navg, 1, 4);
        out.amp_bt.dtype = "uint8";
        out.amp_bt.group = "data_vars";
        out.amp_bt.units = "counts";

        out.prcnt_gd_bt.data = zeros(navg, 1, 4);
        out.prcnt_gd_bt.dtype = "uint8";
        out.prcnt_gd_bt.group = "data_vars";
        out.prcnt_gd_bt.units = "%";

        out.time.data = zeros(navg);
        out.time.dtype = "float64";
        out.time.group = "coords";
        out.time.units = "";

        out.etime_gps.data = zeros(navg);
        out.etime_gps.dtype = "float64";
        out.etime_gps.group = "coords";
        out.etime_gps.units = "";

        out.elatitude_gps.data = zeros(navg);
        out.elatitude_gps.dtype = "float64";
        out.elatitude_gps.group = "data_vars";
        out.elatitude_gps.units = "deg";

        out.elongitude_gps.data = zeros(navg);
        out.elongitude_gps.dtype = "float64";
        out.elongitude_gps.group = "data_vars";
        out.elongitude_gps.units = "deg";

        out.time_gps.data = zeros(navg);
        out.time_gps.dtype = "float64";
        out.time_gps.group = "coords";
        out.time_gps.units = "";

        out.latitude_gps.data = zeros(navg);
        out.latitude_gps.dtype = "float64";
        out.latitude_gps.group = "data_vars";
        out.latitude_gps.units = "deg";

        out.longitude_gps.data = zeros(navg);
        out.longitude_gps.dtype = "float64";
        out.longitude_gps.group = "data_vars";
        out.longitude_gps.units = "deg";

        out.ntime.data = zeros(navg);
        out.ntime.dtype = "float64";
        out.ntime.group = "coords";
        out.ntime.units = "";

        out.flags.data = zeros(navg);
        out.flags.dtype = "float32";
        out.flags.group = "data_vars";
        out.flags.units = "";
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function size = find_filesize()
        fseek(fid, 0, "eof");
        size = ftell(fid);
        frewind(fid);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = init_data(nens)
        out = struct('data_vars', struct(), 'coords', struct(),...
                'attrs', struct(), 'units', struct(), 'sys', struct());
        out.attrs.inst_make = 'TRDI';
        out.attrs.inst_model = 'Workhorse';
        out.attrs.inst_type = 'ADCP';
        out.attrs.rotate_vars = {'vel'};
        % Currently RDI doesn't use IMUs
        out.attrs.has_imu = false;
        fields = fieldnames(ensemble);
        for i = 1:numel(fields)
            key = fields{i};
            group = ensemble.(key).group;
            dtype = ensemble.(key).dtype;
            units = ensemble.(key).units;
            shape = size(ensemble.(key).data);
            shape(1) = nens;
            if startsWith(dtype, "float")
                out.(group).(key) = nan(shape);
            else
                out.(group).(key) = zeros(shape);
            end
            out.units.(key) = units;
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function search_buffer()
        % Check to see if the next bytes indicate the beginning of a
        % data block.  If not, search for the next data block, up to
        % _search_num times.
        id1 = fread(fid,2,'uint8');
        search_cnt = 0;
        while search_cnt < search_num_ && (...
                (id1(1) ~= 127 || id1(2) ~= 127) || ~checkheader())
            search_cnt = search_cnt + 1;
            nextbyte = fread(fid,1,'uint8');
            id1(2) = id1(1);
            id1(1) = nextbyte;
        end
        if search_cnt == search_num_
            ME = MException('MATLAB:read_rdi:search_buffer',...
                'Searched %d entries... Bad data encountered. -> %d,%d',...
                search_cnt, id1(1), id1(2));
            throwAsCaller(ME);
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function skip_Ncol(n)
        fseek(fid, n*cfg.n_cells, 0);
        nbyte_ = 2 + n*cfg.n_cells;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function skip_Nbyte(n)
        fseek(fid, n, 0);
        nbyte_ = 2 + n;
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function check_offset(offset)
        if offset ~= 4 && fixoffset_ == 0
            fixoffset_ = offset - 4;
        end
        fseek(fid, 4 + fixoffset_, 0); 
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = remove_end(data, iens)
        fprintf('  Encountered end of file.  Cleaning up data.\n')
        out = data;
        for qq = 1 : numel(vars_read)
            ky = vars_read{qq};
            try
                group = ensemble.(ky).group;     
                l = length(size((out.(group).(ky))));
                otherdims = repmat({':'},1,l-1);
                out.(group).(ky) = out.(group).(ky)(1:iens-1,otherdims{:});
            catch
                continue
            end            
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = finalize(data)
        out = data;
        % Remove the attributes from the data that were never loaded.
        fields = fieldnames(ensemble);
        for qq = 1 : numel(fields)
            ky = fields{qq};
            if ~any(strcmp(vars_read, ky))   
                try
                    group = ensemble.(ky).group;
                    out.(group) = rmfield(out.(group), ky);                    
                catch
                    continue
                end
                
            end
        end
        fields = fieldnames(cfg);
        for qq = 1 : numel(fields)
            out.attrs.(fields{qq}) = cfg.(fields{qq});
        end
        out.attrs.fs = (data.attrs.sec_between_ping_groups * ...
            data.attrs.pings_per_ensemble)^-1;

        fields = fieldnames(ensemble);
        for qq = 1 : numel(fields)
            ky = fields{qq};
            try
                group = ensemble.(ky).group;
                shape = size(out.(group).(ky));
                if shape(end) == cfg.n_cells
                    out.(group).(ky) = ...
                        permute(out.(group).(ky), [1,2,4,3]);
                end
            catch
                continue
            end            
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function out = remove_gps_duplicates(dat)
        % Removes duplicate and nan timestamp values in 'time_gps' 
        % coordinate, and adds hardware (ADCP DAQ) timestamp corresponding
        % to GPS acquisition (in addition to the GPS unit's timestamp).
        out = dat;
        out.data_vars.hdwtime_gps = out.coords.time;
        out.units.hdwtime = 'seconds since 1970-01-01 00:00:00';

        % Remove duplicate timestamp values, if applicable
        [out.coords.time_gps, idx] = unique(out.coords.time_gps);
     
        % Remove nan values, if applicable
        nans = isnan(out.coords.time_gps);
        out.coords.time_gps = out.coords.time_gps(~nans);

        dv_fields = fieldnames(out.data_vars);
        for qq = 1:numel(dv_fields)
            key = dv_fields{qq};
            if contains(key,'gps')
                out.data_vars.(key) = out.data_vars.(key)(idx);
                if sum(nans) > 0
                    out.data_vars.(key) = out.data_vars.(key)(~nans);
                end
            end
        end
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>    

end

