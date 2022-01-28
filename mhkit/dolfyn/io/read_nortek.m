function ds=read_nortek(filename,options)

%%%%%%%%%%%%%%%%%%%%
%     Read a classic Nortek (AWAC and Vector) datafile 
%     
% Parameters
% ------------
%     filename: string
%         Filename of Nortek file to read.
%     userdata: bool or string (optional)
%         true, false, or string of userdata.json filename (default true)
%         Whether to read the '<base-filename>.userdata.json' file.
%     do_checksum: bool (optional)
%         Whether to perform the checksum of each data block. (default
%         false)
%     nens: nan, int, or 2-element array (optional)
%         nan (default: read entire file), int, or 2-element tuple 
%         (start, stop) Number of pings to read from the file.
%
%     call with options -> read_nortek(filename, userdata=false, do_checksum=true, nens=12) 
%
% Returns
% ---------
%     data: structure 
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
        ME = MException('MATLAB:read_nortek','filename must be a string');
        throw(ME);
    end
    
    % check to see if the file exists
    if ~isfile(filename)
        ME = MException('MATLAB:read_nortek','file does not exist');
        throw(ME);
    end
    
    % check to make sure userdata is bool or string
    if ~isa(options.userdata, 'logical') && ~isa(options.userdata, 'string')
        ME = MException('MATLAB:read_nortek','userdata must be a logical or string');
        throw(ME);
    end
    
    % check to make sure do_checksum is bool
    if ~isa(options.do_checksum, 'logical') 
        ME = MException('MATLAB:read_nortek','do_checksum must be a logical');
        throw(ME);
    end
    
    % check to make sure nens is numeric or nan
    if ~all(isa(options.nens, 'numeric'))
        ME = MException('MATLAB:read_nortek','nens must be numeric or nan');
        throw(ME);
    end
    
    % check to make sure if nens is numeric that its length is equal to 1 or 2
    nstart = 0;
    if ~isnan(options.nens)
        if length(options.nens) < 1 || length(options.nens) > 2
            ME = MException('MATLAB:read_nortek','nens must be a single value or tuple');
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
    endian1 = fread(fid,2,'uint16','b');
    frewind(fid);
    endian2 = fread(fid,2,'uint16','l');
    if all(endian1 == [1445;24])
        endian = 'b';
    elseif all(endian2 == [1445;24])
        endian = 'l';
    else
        % this is not a nortek file. Move to next reader
        ME = MException('MATLAB:read_nortek',['could not determine the ' ...
            'endianness of the file. Are you sure this is a Nortek file?']);
        throw(ME);        
    end
    clearvars endian1 endian2;
    frewind(fid);
    % initialize variables
    config = struct;
    data = struct;
    n_samp_guess = 1;
    filesize = find_filesize();
    c = 1;
    dtypes = {}; 
    lastread = {{},{},{},{},{}};
    vec_data = struct();
    vec_sysdata = struct();
    awac_profile = struct();
    ahrsid_ = struct();
    orient_dnames = {};
    load("nortek_defs.mat", "vec_data", "vec_sysdata", "awac_profile");
    % matlab demands field names begin with a letter so the actual function
    % code is the fieldname without the leading a
    fun_map = struct(  'a0x00',{@read_user_cfg},...
                       'a0x04',{@read_head_cfg},...
                       'a0x05',{@read_hw_cfg},...
                       'a0x07',{@read_vec_checkdata},...
                       'a0x10',{@read_vec_data},...
                       'a0x11',{@read_vec_sysdata},...
                       'a0x12',{@read_vec_hdr},...
                       'a0x71',{@read_microstrain},...
                       'a0x20',{@read_awac_profile},...
                       'sci_awac_profile',{@sci_awac_profile},...
                       'sci_vec_data',{@sci_vec_data},...
                       'sci_vec_sysdata',{@sci_vec_sysdata},...
                       'sci_microstrain',{@sci_microstrain});
    % Read the header:
    if read_id == 5
        read_hw_cfg;
    else
        ME = MException('MATLAB:read_nortek',['Corrupted block found in ' ...
            'header']);
        throw(ME);
    end
    if read_id == 4
        read_head_cfg;
    else
        ME = MException('MATLAB:read_nortek',['I/O error: The file does ' ...
            'not appear to be a Nortek data file.']);
        throw(ME);
    end
    if read_id == 0
        read_user_cfg;
    else
        ME = MException('MATLAB:read_nortek',['I/O error: The file does ' ...
            'not appear to be a Nortek data file.']);
        throw(ME);
    end
    serial_num = char(config.('serialNum'));
    serial_num = upper(serial_num(1:3));
    pnow = ftell(fid);
    % Run the appropriate initialization routine (e.g. init_ADV).
    if serial_num == "WPR"
        config.('config_type') = "AWAC";
        init_AWAC();
    elseif serial_num == "VEC"
        config.('config_type') = "ADV";
        init_ADV();
    end
    if ~isnan(npings)
        n_samp_guess = npings + 1;
    end
    % return the file position to pnow
    fseek(fid,pnow,"bof");
    % finalize initialize data(attrs)
    if config.('NBurst') > 0
        data.('attrs').("DutyCycle_NBurst") = config.('NBurst');
        data.('attrs').("DutyCycle_NCycle") = config.('MeasInterval') * ...
            config.('fs');
    end
    burst_start = false([n_samp_guess,1]);
    data.('attrs').('fs') = config.('fs');
    if config.('coord_sys_axes') == 'XYZ'        
        data.('attrs').('coord_sys') = 'inst';
    elseif config.('coord_sys_axes') == 'ENU'        
        data.('attrs').('coord_sys') = 'earth';
    elseif config.('coord_sys_axes') == 'beam'        
        data.('attrs').('coord_sys') = 'beam';
    end
    data.('attrs').("has_imu") = 0;
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        End Initialize
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Read File
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    readfile();
    fclose(fid);    % close the file 

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                       Post Process
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    dat2sci();

    rotmat = nan;
    declin = nan;
    fn = fieldnames(userdata);
    for k=1:numel(fn)
        if( contains(fn{k},'rotmat') )
            rotmat = userdata.(fn{k});        
        elseif( contains(fn{k},'dec') )
            declin = userdata.(fn{k});
        else
            data.attrs.(fn{k}) = userdata.(fn{k});
        end
    end

    handle_nan();
    ds = create_dataset(data);
    % This next call is very specific to xarray so it is being skipped
    % ds = set_coords(ds, ds.coord_sys);
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
        ds.orientmat.dims = {'earth', 'inst', 'time'};
        ds.orientmat.coords.earth = {'E' 'N' 'U'};
        ds.orientmat.coords.inst = {'X' 'Y' 'Z'};
        ds.orientmat.coords.time = ds.time;
    end

    if ~isnan(rotmat)
        ds = set_inst2head_rotmat(ds, rotmat);
    end
    if ~isnan(declin)
        ds = set_declination(ds, declin);
    end
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        End Post Process
    %                        End Function
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Reader Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    function userdata = read_userdata(filename, userdata)
        %%%%%%%%%%%%%%%%%%%%
        %     Reads a userdata.json file and returns the data it contains
        %     as a structure
        %     
        % Parameters
        % ------------
        %     filename: string
        %         Filename of Nortek file to read.
        %     userdata: bool or string 
        %         true, false, or string of userdata.json filename         
        %
        % Returns
        % ---------
        %     userdata: structure          
        %        
        %%%%%%%%%%%%%%%%%%%%
        if isa(userdata, 'logical')
            if ~userdata
                userdata = struct;
                return
            else
                % determine name of userdata file from base name
                basefile = split(filename,".");
                jsonfile = basefile(1) + ".userdata.json";
            end
        else
            % use name directly as supplied
            jsonfile = userdata;
        end
        % make sure the file exists
        if ~isfile(jsonfile)
            userdata = struct;
            return
        end
        % read the json data
        fid = fopen(jsonfile); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        userdata = jsondecode(str);

        % quality checks for userdata

        % if the following fields exist, rename them
        nm_list = {'body2head_rotmat', 'body2head_vec'};
        for i = 1 : length(nm_list)
            nm = nm_list{i};
            if isfield(userdata, nm)
                new_name = append('inst', nm(5:end));
                userdata.(new_name) = userdata.(nm);                
                userdata = rmfield(userdata, nm);
            end
        end
        % if inst2head_rotmat data = identity, eye, 1, or 1. then change
        % data to be a 3x3 identity matrix
        if isfield(userdata, 'inst2head_rotmat')
            % check if its a string
            if isa(userdata.('inst2head_rotmat'), 'char')
                if strcmp(userdata.('inst2head_rotmat'),'identity') ||...
                    strcmp(userdata.('inst2head_rotmat'),'eye')
                    userdata.('inst2head_rotmat') = eye(3);
                end
            else % if not maybe its numeric check if its length 1
                if userdata.('inst2head_rotmat')
                    if length(userdata.('inst2head_rotmat')) == 1
                        if userdata.('inst2head_rotmat') == 1
                            userdata.('inst2head_rotmat') = eye(3);
                        end
                    end
                end
            end            
        end
        % Make sure that coord_sys is not in the userdata
        if isfield(userdata, 'coord_sys')
            msgtext = ['The instrument coordinate system (coord_sys) should' ...
                ' not be specified in the .userdata.json file, remove this' ...
                ' and read the file again.'];
            ME = MException('MATLAB:read_nortek:read_userdata',msgtext);
            throwAsCaller(ME)
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function readfile()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Loops through the file to read data              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        retval = false;
        nlines = nan;
        try
            while ~retval
                if c == nlines
                    break
                end
                retval = readnext();
                if retval == 10
                    findnext(true);
                    retval = false;
                end
                if ~isnan(npings) && c >= nprings
                    if any(strcmp(dtypes,'microstrain'))
                        try
                            readnext();
                        catch
                            continue
                        end
                    end
                break
                end
            end
        catch 
            fprintf("End of file at %d bytes.\n",ftell(fid))
        end
        c = c - 1;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = readnext()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Reads the next id in the file              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        id = read_id();
        id = dec2hex(id);
        if length(id)==1
            id = append("0" + id);
        end
        id = append("a0x" + id);
        if isfield(fun_map,id)
            out = feval(fun_map.(id));
            fun_name = extractAfter(func2str(fun_map.(id)),...
                "read_nortek/read_");
            lastread =  [{fun_name},lastread(1:end-1)];
        else
            warning(append('Unrecognized identifier: ', ...
                extractAfter(id,"a")))
            fseek(fid,-2,0);
            out = 10;
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function id = read_id()
        %%%%%%%%%%%%%%%%%%%%
        %     Read the next 'ID' from the file.              
        %
        % Returns
        % ---------
        %     id: int          
        %        
        %%%%%%%%%%%%%%%%%%%%
        read_bytes = fread(fid,2,'uchar',endian);
        if read_bytes(1) ~= 165
            % corrupted block            
%             warning("Corrupted data block sync code (%d, %d) found" + ...
%                 " in ping %d. Searching for next valid code", read_bytes(1),...
%                 read_bytes(2), c);
            id = findnext(false);
            fseek(fid, 2, 0);            
        else
            id = read_bytes(2);
        end        
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function read_hw_cfg()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read hardware configuration              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        fseek(fid, 2, 0); % skip size of structure
        temp = fread(fid,8,'char',endian);
        config.('serialNum') = convertCharsToStrings(char(temp));
        config.('ProLogID') = fread(fid,1,'int16',endian);
        temp = fread(fid,4,'char',endian);
        config.('ProLogFWver') = convertCharsToStrings(char(temp));
        temp = fread(fid,6,'uint16',endian);
        config.('config') = temp(1);
        config.('freq') = temp(2);
        config.('PICversion') = temp(3);
        config.('HWrevision') = temp(4);
        config.('recSize') = temp(5) * 65536;
        config.('status') = temp(6);
        fseek(fid, 12, 0); % skip spare
        config.('FWversion') = fread(fid,1,'uint32',endian);
        do_checksum;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function read_head_cfg()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read head configuration              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        fseek(fid, 2, 0);   % skip size of structure
        fseek(fid, 2, 0);   % skip config
        config.('freq') = fread(fid,1,'ushort',endian);
        fseek(fid, 2, 0);   % skip head type
        fseek(fid, 12, 0);  % skip head serial number
        fseek(fid, 8, 0);   % skip beginning of system data
        temp = fread(fid,9,'int16',endian);
        temp = temp./4096.;
        config.('beam2inst_orientmat') = transpose(reshape(temp,[3,3]));
        fseek(fid, 150, 0); % skip the rest of system data
        fseek(fid, 22, 0);  % skip spare
        fseek(fid, 2, 0);   % skip number of beams
        do_checksum;
    end 

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function read_user_cfg()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read user configuration             
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % 
        Profile_Timing = {'single', 'continuous'};
        sync_out = {'middle', 'end'};
        coord_sys_axes = {'ENU', 'XYZ', 'beam'};
        output_format = {'Vector', 'ADV'};
        vel_scale = [1, 0.1];
        mode_test_filter_output = {'total', 'correction_only'};
        rate = {'1hz', '2hz'};
        cell_position = {'fixed', 'dynamic'};
        dynamic_pos_type = {'pct of mean press', 'pct of min re'};
        %
        fseek(fid, 2, 0);   % skip size of structure
        config.('Transmit').('pulse_length') = ...
            fread(fid,1,'ushort',endian);
        config.('Transmit').('blank_distance') = ...
            fread(fid,1,'ushort',endian);
        config.('Transmit').('receive_length') = ...
            fread(fid,1,'ushort',endian);
        config.('Transmit').('time_between_pings') = ...
            fread(fid,1,'ushort',endian);
        config.('Transmit').('time_between_bursts') = ...
            fread(fid,1,'ushort',endian);
        config.('Npings') = fread(fid,1,'ushort',endian);
        config.('AvgInterval') = fread(fid,1,'ushort',endian);  
        config.('NBeams') = fread(fid,1,'ushort',endian); 
        config.('TimCtrlReg') = int2binarray(...
            fread(fid,1,'ushort',endian),16); 
        % From the nortek system integrator manual
        % (note: bit numbering is zero-based)
        treg = int16(config.('TimCtrlReg'));
        config.('Profile_Timing') = Profile_Timing{treg(2)+1};
        config.('Burst_Mode') = logical(~treg(3));
        config.('Power_Level') = treg(6) + 2 * treg(7) + 1;
        config.('sync_out') = sync_out{treg(8)+1};
        config.('Sample_on_Sync') = logical(treg(9));
        config.('Start_on_Sync') = logical(treg(10));
        config.('PwrCtrlReg') = int2binarray(...
            fread(fid,1,'ushort',endian),16); 
        config.('A1') = fread(fid,1,'ushort',endian); 
        config.('B0') = fread(fid,1,'ushort',endian); 
        config.('B1') = fread(fid,1,'ushort',endian); 
        config.('CompassUpdRate') = fread(fid,1,'ushort',endian); 
        config.('coord_sys_axes') = ...
            coord_sys_axes{fread(fid,1,'ushort',endian) + 1}; 
        config.('NBins') = fread(fid,1,'ushort',endian); 
        config.('BinLength') = fread(fid,1,'ushort',endian); 
        config.('MeasInterval') = fread(fid,1,'ushort',endian);
        temp = fread(fid,6,'char',endian);
        config.('DeployName') = convertCharsToStrings(char(temp(1:4)));
        config.('WrapMode') = fread(fid,1,'ushort',endian);
        config.('ClockDeploy') = fread(fid,3,'ushort',endian);
        config.('DiagInterval') = fread(fid,1,'ulong',endian);
        config.('Mode0') = int2binarray(fread(fid,1,'ushort',endian), 16);
        config.('AdjSoundSpeed') = fread(fid,1,'ushort',endian); 
        config.('NSampDiag') = fread(fid,1,'ushort',endian); 
        config.('NBeamsCellDiag') = fread(fid,1,'ushort',endian); 
        config.('NPingsDiag') = fread(fid,1,'ushort',endian); 
        config.('ModeTest') = ...
            int2binarray(fread(fid,1,'ushort',endian), 16); 
        config.('AnaInAddr') = fread(fid,1,'ushort',endian); 
        config.('SWVersion') = fread(fid,1,'ushort',endian); 
        fseek(fid, 2, 0);   % skip salinity
        config.('VelAdjTable') = fread(fid,90,'ushort',endian);
        temp = fread(fid,80,'char',endian);
        config.('Comments') = convertCharsToStrings(char(temp)); 
        fseek(fid, 100, 0);   % skip spare, prcoessing method, spare
        config.('Mode1') = int2binarray(fread(fid,1,'ushort',endian), 16);
        config.('DynPercPos') = fread(fid,1,'ushort',endian);
        config.('T1w') = fread(fid,1,'ushort',endian);
        config.('T2w') = fread(fid,1,'ushort',endian);
        config.('T3w') = fread(fid,1,'ushort',endian);
        config.('NSamp') = fread(fid,1,'ushort',endian);
        fseek(fid, 4, 0);   % skip A1 & B0
        config.('NBurst') = fread(fid,1,'ushort',endian);
        fseek(fid, 2, 0);   % skip spare
        config.('AnaOutScale') = fread(fid,1,'ushort',endian);
        config.('CorrThresh') = fread(fid,1,'ushort',endian);
        fseek(fid, 2, 0);   % skip spare
        config.('TiLag2') = fread(fid,1,'ushort',endian);
        fseek(fid, 30, 0);   % skip spare
        config.('QualConst')  = fread(fid,8,'ushort',endian);
        do_checksum;
        config.('mode').('user_sound') = config.('Mode0')(1);
        config.('mode').('diagnostics_mode') = config.('Mode0')(2);
        config.('mode').('analog_output_mode') = config.('Mode0')(3);
        config.('mode').('output_format') = ...
            output_format{config.('Mode0')(4)+1};
        config.('mode').('vel_scale') = vel_scale(config.('Mode0')(5)+1);
        config.('mode').('serial_output') = config.('Mode0')(6);
        config.('mode').('reserved_EasyQ') = config.('Mode0')(7);
        config.('mode').('stage') = config.('Mode0')(8);
        config.('mode').('output_power') = config.('Mode0')(9);
        config.('mode').('mode_test_use_DSP') = config.('ModeTest')(1);
        config.('mode').('mode_test_filter_output') = ...
            mode_test_filter_output{config.('ModeTest')(2) +1}; % noqa
        config.('mode').('rate') = rate{config.('Mode1')(1)+1};
        config.('mode').('cell_position') = cell_position{...
            config.('Mode1')(2)+1};
        config.('mode').('dynamic_pos_type') = dynamic_pos_type{...
            config.('Mode1')(3)+1}; % noqa
    end 

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_vec_checkdata()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x07              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        checknow = struct();
        fseek(fid, 2, 0);   % skip size of structure
        checknow.Samples = fread(fid,1,'ushort',endian);
        n = checknow.Samples;
        checknow.First_samp = fread(fid,1,'ushort',endian);
        checknow.Amp1 = fread(fid,n,'uint8',endian);
        checknow.Amp2 = fread(fid,n,'uint8',endian);
        checknow.Amp3 = fread(fid,n,'uint8',endian);
        do_checksum()        
        if ~isfield(data.attrs.config,'checkdata')
            data.attrs.config.checkdata = checknow;
        else
            temp = concat_struct(data.attrs.config.checkdata, checknow);
            data.attrs.config.checkdata = temp;
        end 
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_vec_data()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x10              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if ~isfield(data.data_vars,'vel')
            init_data(vec_data);
            if ~any(strcmp(dtypes,'vec_data'))
                dtypes{end+1} = 'vec_data';
            end             
        end

        data.sys.AnaIn2LSB(c) = fread(fid,1,'uint8',endian);
        data.sys.Count(c) = fread(fid,1,'uint8',endian);
        data.data_vars.PressureMSB(c) = fread(fid,1,'uint8',endian);
        data.sys.AnaIn2MSB(c) = fread(fid,1,'uint8',endian);
        data.data_vars.PressureLSW(c) = fread(fid,1,'ushort',endian);
        data.sys.AnaIn1(c) = fread(fid,1,'ushort',endian);
        data.data_vars.vel(c,1,1) = fread(fid,1,'short',endian);
        data.data_vars.vel(c,1,2) = fread(fid,1,'short',endian);
        data.data_vars.vel(c,1,3) = fread(fid,1,'short',endian);
        data.data_vars.amp(c,1,1) = fread(fid,1,'uint8',endian);
        data.data_vars.amp(c,1,2) = fread(fid,1,'uint8',endian);
        data.data_vars.amp(c,1,3) = fread(fid,1,'uint8',endian);
        data.data_vars.corr(c,1,1) = fread(fid,1,'uint8',endian);
        data.data_vars.corr(c,1,2) = fread(fid,1,'uint8',endian);
        data.data_vars.corr(c,1,3) = fread(fid,1,'uint8',endian);

        do_checksum;
        c = c + 1;
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_vec_sysdata()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x11             
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if lastread{1} == "vec_checkdata" && lastread{2} == "vec_hdr"
            burst_start(c) = true;
        end
        if ~isfield(data.coords,'time')
            init_data(vec_sysdata);
            if ~any(strcmp(dtypes,'vec_sysdata'))
                dtypes{end+1} = 'vec_sysdata';
            end             
        end
        fseek(fid, 2, 0); % skip size of structure    
        data.coords.time(c) = rd_time();
        data.data_vars.batt(c) = fread(fid,1,'ushort',endian);
        data.data_vars.c_sound(c) = fread(fid,1,'ushort',endian);
        data.data_vars.heading(c) = fread(fid,1,'short',endian);        
        data.data_vars.pitch(c) = fread(fid,1,'short',endian);        
        data.data_vars.roll(c) = fread(fid,1,'short',endian);
        data.data_vars.temp(c) = fread(fid,1,'ushort',endian);
        data.data_vars.error(c) = fread(fid,1,'uint8',endian);
        data.data_vars.status(c) = fread(fid,1,'uint8',endian);
        data.sys.AnaIn(c) = fread(fid,1,'ushort',endian);
        do_checksum();
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_vec_hdr()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x12              
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        hdrnow = struct();
        fseek(fid, 2, 0); % skip size of structure        
        hdrnow.time = rd_time();
        hdrnow.NRecords = fread(fid,1,'ushort',endian);
        hdrnow.Noise1   = fread(fid,1,'uint8',endian);
        hdrnow.Noise2   = fread(fid,1,'uint8',endian);
        hdrnow.Noise3   = fread(fid,1,'uint8',endian); 
        hdrnow.Spare0   = fread(fid,1,'uint8',endian);
        hdrnow.Corr1    = fread(fid,1,'uint8',endian);
        hdrnow.Corr2    = fread(fid,1,'uint8',endian);
        hdrnow.Corr3    = fread(fid,1,'uint8',endian);
        hdrnow.Spare1   = fread(fid,21,'char',endian);
        do_checksum();
        if ~isfield(data.attrs.config,'data_header')
            data.attrs.config.data_header = hdrnow;
        else
            temp = concat_struct(data.attrs.config.data_header, hdrnow);
            data.attrs.config.data_header = temp;
        end
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_microstrain()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x71             
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if c == 1
            warning(['First "microstrain data" block is before first ' ...
                '"vector system data" block.'])
        else
            c = c - 1;
        end
        fseek(fid, 3, 0);
        ahrsid = fread(fid,1,'uint8',endian);
        if isfield(ahrsid_, 'ahrsid') && ahrsid ~= ahrsid_.ahrsid
            warning('AHRS_ID changes mid-file!')
        end

        if any([195, 204, 210, 211] == ahrsid)
            ahrsid_.ahrsid = ahrsid;
        end

        data.attrs.has_imu = 1;

        if ~isfield(data.data_vars,'accel')
            dtypes{end+1} = 'microstrain';
            switch ahrsid
                case 195
                    orient_dnames = {"accel", "angrt", "orientmat"};
                    data.data_vars.accel = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.angrt = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.orientmat = nan([n_samp_guess,1,3,3],...
                        'single');
                    rv = {'accel', 'angrt'};
                    for i = 1:length(rv)
                        if ~any(strcmp(data.attrs.rotate_vars,rv{i}))
                            data.attrs.rotate_vars{end+1} = rv{i};
                        end 
                    end
                    data.units.accel = "m/s^2";
                    data.units.angrt = "rad/s";
                case {204, 210}
                    orient_dnames = {"accel", "angrt", "mag", "orientmat"};
                    data.data_vars.accel = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.angrt = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.mag = nan([n_samp_guess,1,3],...
                        'single');
                    rv = {'accel', 'angrt', 'mag'};
                    for i = 1:length(rv)
                        if ~any(strcmp(data.attrs.rotate_vars,rv{i}))
                            data.attrs.rotate_vars{end+1} = rv{i};
                        end 
                    end
                    if ahrsid == 204
                        data.data_vars.orientmat = nan([n_samp_guess,1,3,3],...
                        'single');
                    end
                    data.units.accel = "m/s^2";
                    data.units.angrt = "rad/s";
                    data.units.mag   = "gauss";
                case 211
                    orient_dnames = {"angrt", "accel", "angrt"};
                    data.data_vars.accel = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.angrt = nan([n_samp_guess,1,3],...
                        'single');
                    data.data_vars.mag = nan([n_samp_guess,1,3],...
                        'single');
                    rv = {'accel', 'angrt', 'mag'};
                    for i = 1:length(rv)
                        if ~any(strcmp(data.attrs.rotate_vars,rv{i}))
                            data.attrs.rotate_vars{end+1} = rv{i};
                        end 
                    end
                    data.units.accel = "m/s^2";
                    data.units.angrt = "rad/s";
                    data.units.mag   = "gauss";
                otherwise
            end
        end

        switch ahrsid
            case 195
                temp = fread(fid,15,'single',endian);
                data.data_vars.angrt(c,:,:) = temp(1:3);
                data.data_vars.accel(c,:,:) = temp(4:6);
                data.data_vars.orientmat(c,:,:,1) = temp( 7: 9);
                data.data_vars.orientmat(c,:,:,2) = temp(10:12);
                data.data_vars.orientmat(c,:,:,3) = temp(13:15);
                fseek(fid, 4, 0);
            case 204
                temp = fread(fid,18,'single',endian);
                data.data_vars.accel(c,:,:) = temp(1:3);
                data.data_vars.angrt(c,:,:) = temp(4:6);
                data.data_vars.mag(c,:,:)   = temp(7:9);
                data.data_vars.orientmat(c,:,:,1) = temp(10:12);
                data.data_vars.orientmat(c,:,:,2) = temp(13:15);
                data.data_vars.orientmat(c,:,:,3) = temp(16:18);
                fseek(fid, 6, 0);
            case 211
                temp = fread(fid,9,'single',endian);
                data.data_vars.angrt(c,:,:) = temp(1:3);
                data.data_vars.accel(c,:,:) = temp(4:6);
                data.data_vars.mag(c,:,:)   = temp(7:9);
                fseek(fid, 6, 0);
            otherwise
                warning("Unrecognized IMU identifier: %d", ahrsid);
                fseek(fid,-2,0);
                out = 10;
        end
        do_checksum();
        c = c + 1;
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = read_awac_profile()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read Code: 0x20 (AWAC velocity data)            
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        nbins = config.NBins;
        if ~isfield(data.data_vars,'temp')
            init_data(awac_profile);
            if ~any(strcmp(dtypes,'awac_profile'))
                dtypes{end+1} = 'awac_profile';
            end 
        end
        % Note: docs state there is 'fill' byte at the end, if nbins is odd,
        % but doesn't appear to be the case
        n = config.NBeams;
        fseek(fid, 2, 0); % skip size of structure        
        data.coords.time(c) = rd_time();
        data.data_vars.error(c) = fread(fid,1,'ushort',endian);
        data.sys.AnaIn1(c) = fread(fid,1,'ushort',endian);
        data.data_vars.batt(c) = fread(fid,1,'ushort',endian);
        data.data_vars.c_sound(c) = fread(fid,1,'ushort',endian);
        data.data_vars.heading(c) = fread(fid,1,'ushort',endian);
        data.data_vars.pitch(c) = fread(fid,1,'ushort',endian);
        data.data_vars.roll(c) = fread(fid,1,'ushort',endian);
        p_msb = fread(fid,1,'char*1',endian);
        data.data_vars.status(c) = fread(fid,1,'*char',endian);
        p_lsw = fread(fid,1,'ushort',endian);
        data.data_vars.temp(c) = fread(fid,1,'ushort',endian);
        data.data_vars.pressure(c) = (65536 * p_msb + p_lsw);
        % The nortek system integrator manual specifies an 88byte 'spare'
        fseek(fid, 88, 0); % skip size of structure 
        temp_vel = fread(fid,n * nbins,'int16',endian);
        temp_amp = fread(fid,n * nbins,'uint8',endian);
        for i = 1:n
            data.data_vars.vel(c,1,i,:) = temp_vel((i-1)*nbins+1:i*nbins);
            data.data_vars.amp(c,1,i,:) = temp_amp((i-1)*nbins+1:i*nbins);
        end
        do_checksum();
        c = c + 1;
        out = false;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    %% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %                        Helper Functions
    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function init_ADV()
        config.('fs') = 512 / config.('AvgInterval');
        data.('attrs').('config') = config;
        data.('attrs').('inst_make') = 'Nortek';
        data.('attrs').('inst_model') = 'Vector';
        data.('attrs').('inst_type') = 'ADV';
        data.('attrs').('rotate_vars') = {'vel'};
        data.('attrs').('freq') = config.('freq');
        data.('attrs').('SerialNum') = config.('serialNum');
        data.('attrs').('Comments') = config.('Comments');
        data.('data_vars').('beam2inst_orientmat') = ...
            config.('beam2inst_orientmat');
        data.('coords') = struct;
        data.('units')  = struct;
        data.('sys')    = struct;
        dlta = code_spacing('0x11');
        n_samp_guess = int32(filesize / dlta + 1);
        n_samp_guess = n_samp_guess * int32(config.("fs"));
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function init_AWAC()
        config.('fs') = 1. / config.('AvgInterval');
        data.('attrs').('config') = config;
        data.('attrs').('inst_make') = 'Nortek';
        data.('attrs').('inst_model') = 'AWAC';
        data.('attrs').('inst_type') = 'ADCP';
        data.('attrs').('rotate_vars') = {'vel'};
        data.('attrs').('freq') = config.('freq');
        data.('attrs').('SerialNum') = config.('serialNum');
        data.('attrs').('Comments') = config.('Comments');
        data.('attrs').('n_beams') = config.('NBeams');
        data.('attrs').('avg_interval') = config.('AvgInterval');
        data.('data_vars').('beam2inst_orientmat') = ...
            config.('beam2inst_orientmat');
        data.('coords') = struct;
        data.('units')  = struct;
        data.('sys')    = struct;
        space = code_spacing('0x20');
        if space == 0
            n_samp_guess = 1;
        else
            n_samp_guess = int32(filesize / space + 1);
        end
    end 

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function do_checksum()
        if options.do_checksum
            warning("Not implemented")
        else
            fseek(fid,2,0);
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function loc = findnextid(id)
        if isa(id,'char')
            id = hex2dec(id);
        end
        nowid = -1;
        while nowid ~= id
            nowid = read_id();
            if nowid == 16
                shift = 22;
            else
                sz = 2 * fread(fid,1,'ushort',endian);
                shift = sz - 4;
            end
            fseek(fid,shift,0);
        end
        loc = ftell(fid);
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function id = findnext(do_cs)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Find the next data block by checking
        % the checksum and the sync byte(0xa5)          
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        sum = hex2dec('0xb58c');
        cs = 0;
        while ftell(fid) ~= filesize
            val = fread(fid,1,'ushort',endian);
            fseek(fid,-2,0);
            read_bytes = fread(fid,2,'uchar',endian);
            if (all(read_bytes(1) == 165)) && (~do_cs || cs==sum)
                id = read_bytes(2);
                fseek(fid,-2,0);
                break
            end
            sum = sum + cs;
            cs = val;
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function spacing = code_spacing(searchcode)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Find the spacing, in bytes, between a specific hardware code.
        % Repeat this * iternum * times(default 50).
        % Returns the average spacing, in bytes, between the code.            
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        p0 = findnextid(searchcode);
        for i = 1:51
            try
                findnextid(searchcode);
            catch
                break
            end
        end
        spacing = (ftell(fid) - p0) / (i);
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = int2binarray(val, n)
        %%%%%%%%%%%%%%%%%%%%
        %        
        % Parameters
        % ------------
        %     val: integer
        %         integer to convert.
        %     n: integer
        %         number of bits         
        %
        % Returns
        % ---------
        %     out: array of boolean           
        %        
        %%%%%%%%%%%%%%%%%%%%
        out = false(1,n);
        for i = 1 : n
            out(i) = bitand(val,2^(i-1));
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function size = find_filesize()
        fseek(fid, 0, "eof");
        size = ftell(fid);
        frewind(fid);
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function init_data(varstruct)
        %%%%%%%%%%%%%%%%%%%%
        %        
        % Parameters
        % ------------
        %     varstruct: structure
        %         The variable definitions loaded in the init phase that               
        %         describe how to initialize each data variable.
        %         options ("vec_data", "vec_sysdata", "awac_profile")
        %%%%%%%%%%%%%%%%%%%%
        shape_args.("n") = n_samp_guess;
        try
            shape_args.("nbins") = config.("NBins");
        catch
            % pass
        end
        fn = fieldnames(varstruct);
        for i=1:numel(fn)
            nm = fn{i};
            va = varstruct.(nm);
            if va.("group") == ""
                if ~isfield(data,nm)
                    data.(nm) =...
                        empty_array(shape_args, va.("dims"), va.("dtype"));
                end
            else
                if ~isfield(data.(va.("group")),nm)
                    data.(va.("group")).(nm) = ...
                        empty_array(shape_args, va.("dims"), va.("dtype"));
                    data.("units").(nm) = va.("units");
                end
            end
        end       
        
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function dat2sci()
        for i = 1:length(dtypes)
            nm = dtypes{i};
            id = append("sci_" + nm);
            feval(fun_map.(id));
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function sci_awac_profile()
        sci_data(awac_profile);
        % Calculate the ranges.
        if config.freq == 2000
            cs_coef = 0.0239;
        elseif config.freq == 1000
            cs_coef = 0.0478;
        elseif config.freq == 600
            cs_coef = 0.0797;
        else
            cs_coef = 0.1195;
        end
        h_ang = 25 * (pi / 180);  % Head angle is 25 degrees for all awacs.
        % Cell size
        cs = round(single(config.BinLength)/256. * cs_coef * cos(h_ang),2);
        % Blanking distance
        bd = round(config.Transmit.blank_distance*0.0229*cos(h_ang)-cs,2);

        r = single([1:config.NBins])*cs + bd;
        data.coords.range = r;
        data.attrs.cell_size = cs;
        data.attrs.blank_dist = bd;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function sci_vec_data()
        sci_data(vec_data);
        data.data_vars.pressure = (single(data.data_vars.PressureMSB)*...
            65536 + single(data.data_vars.PressureLSW) ) / 1000.;
        data.units.pressure = "dbar";
        
        data.data_vars = rmfield(data.data_vars, "PressureMSB");
        data.data_vars = rmfield(data.data_vars, "PressureLSW");

        data.data_vars.vel = data.data_vars.vel*config.mode.vel_scale;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function sci_vec_sysdata()
        % Translate the data in the vec_sysdata structure into 
        % scientific units.
        sci_data(vec_sysdata);
        data.sys.sysi = ~isnan(data.coords.time);
        nburst = config.NBurst;
        fs = data.attrs.fs;
        data.data_vars.orientation_down = ~false(size(data.coords.time));
        if nburst == 0
            num_bursts = 1;
            nburst = length(data.coords.time);
        else
            num_bursts = int32(floor(length(data.coords.time)/nburst));
        end
        for nb = 1:(num_bursts)
            iburst = [(nb-1) * nburst + 1; nb * nburst];
            sysi = data.sys.sysi(iburst(1):iburst(2));
            if isempty(sysi)
                break
            end
            inds = find(sysi);
            inds = inds(2:end);
            time_slice = data.coords.time(iburst(1):iburst(2));
            arng = linspace(1,length(time_slice),length(time_slice));
            if length(inds) >= 2
                p = polyfit(inds,time_slice(inds),1);
                data.coords.time(iburst(1):iburst(2)) = polyval(p,arng);
            elseif length(inds) == 1
                data.coords.time(iburst(1):iburst(2)) = ((arng - inds(1))...
                    / (fs * 3600 * 24) + time_slice(inds(1)));
            else
                data.coords.time(iburst(1):iburst(2)) = ...
                    (time_slice(1) + arng/(fs * 3600 * 24));
            end
            tmpd = nan(size(data.data_vars.heading(iburst(1):...
                iburst(2))),"single");
            % The first status bit should be the orientation.
            stat_slc = data.data_vars.status(iburst(1):iburst(2));
            tmpd(sysi) = bitand(stat_slc(sysi),1);
            tmpd = fillgaps(tmpd, true);
            tmpd(isnan(tmpd)) = 0;
            slope = diff(tmpd);
            tmpd(slope<0) = 1;
            tmpd(slope>0) = 0;
            data.data_vars.orientation_down(iburst(1):iburst(2)) = ...
                logical(tmpd);
        end
        data.data_vars.batt = interpgaps(data.data_vars.batt,...
            data.coords.time, false);
        data.data_vars.c_sound = interpgaps(data.data_vars.c_sound,...
            data.coords.time, false);
        data.data_vars.heading = interpgaps(data.data_vars.heading,...
            data.coords.time, false);
        data.data_vars.pitch = interpgaps(data.data_vars.pitch,...
            data.coords.time, false);
        data.data_vars.roll = interpgaps(data.data_vars.roll,...
            data.coords.time, false);
        data.data_vars.temp = interpgaps(data.data_vars.temp,...
            data.coords.time, false);
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function sci_microstrain()
        ormat = false;
        for i = 1:length(orient_dnames)
            % Rotate the MS orientation data (in MS coordinate system)
            % to be consistent with the ADV coordinate system.
            % (x,y,-z)_ms = (z,y,x)_adv
            nm = orient_dnames{i};
            if nm == "orientmat"
                ormat = true;
                continue
            end
            temp = zeros(size(data.data_vars.(nm)));
            temp(:,:,1) = data.data_vars.(nm)(:,:,3) * -1.;
            temp(:,:,2) = data.data_vars.(nm)(:,:,2);
            temp(:,:,3) = data.data_vars.(nm)(:,:,1);
            data.data_vars.(nm) = temp;
        end
        if ormat
            temp = zeros(size(data.data_vars.orientmat));
            temp(:,:,:,1) = data.data_vars.orientmat(:,:,:,3) * -1.;
            temp(:,:,:,2) = data.data_vars.orientmat(:,:,:,2);
            temp(:,:,:,3) = data.data_vars.orientmat(:,:,:,1);
            data.data_vars.orientmat = temp;
            % MS coordinate system is in North-East-Down (NED),
            % we want East-North-Up (ENU)
            temp = zeros(size(data.data_vars.orientmat));
            temp(:,:,1,:) = data.data_vars.orientmat(:,:,2,:);
            temp(:,:,2,:) = data.data_vars.orientmat(:,:,1,:);
            temp(:,:,3,:) = data.data_vars.orientmat(:,:,3,:) * -1;  
            data.data_vars.orientmat = temp;
        end
        if isfield(data.data_vars,"accel")
            data.data_vars.accel = data.data_vars.accel * 9.80665;
        end
        if ahrsid_.ahrsid == 195 || ahrsid_.ahrsid == 211
            data.data_vars.angrt = data.data_vars.angrt * config.fs;
            data.data_vars.accel = data.data_vars.accel * config.fs;
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function sci_data(varstruct)
        % Convert the data to scientific units accordint to vardict.
        % 
        % Parameters
        % ----------
        %     varstruct: structure
        %         The variable definitions loaded in the init phase that               
        %         describe how to initialize each data variable.
        %         options ("vec_data", "vec_sysdata", "awac_profile")
        fn = fieldnames(varstruct);
        for i=1:numel(fn)
            nm = fn{i};
            va = varstruct.(nm);
            if va.("group") == ""
                if isfield(va,'offset')
                    data.nm = data.nm + va.offset;
                end
                if isfield(va,'factor')
                    data.nm = data.nm * va.factor;
                end
            else
                if isfield(va,'offset')
                    data.(va.group).(nm) = data.(va.group).(nm) +...
                        va.offset;
                end
                if isfield(va,'factor')
                    data.(va.group).(nm) = data.(va.group).(nm)*va.factor;
                end
            end            
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    function out = empty_array(shape_args, dims, dtype)
        if isempty(dims)
            shape = [shape_args.("n"), 1];
        else
            shape = [shape_args.("n"), 1, dims{1}];
            if any(strcmp(dims,"nbins")) 
                shape(4) = shape_args.("nbins");
            end
        end
        if dtype == "float32" || dtype == "float64"
            if dtype == "float32"
                dtype = "single";
            else
                dtype = "double";
            end
            out = NaN(shape,dtype);
        else
            out = zeros(shape,dtype);
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function time = rd_time()
        minutes = bcd2char(fread(fid,1,'char*1',endian));
        seconds = bcd2char(fread(fid,1,'char*1',endian));
        day     = bcd2char(fread(fid,1,'char*1',endian));
        hour    = bcd2char(fread(fid,1,'char*1',endian)); 
        year    = bcd2char(fread(fid,1,'char*1',endian)); 
        month   = bcd2char(fread(fid,1,'char*1',endian));   
        if year < 100
            year = year + 1900 + 100 * (year < 90);
        end
        format long
        time = datetime(year,month,day,hour,minutes,seconds);
        time = datenum(datestr(time), 'dd-mmm-yyyy HH:MM:SS');
        % to convert back to datetime
        % date = datetime(datenum value,'ConvertFrom','datenum')
        % date = date.Format = 'dd-MMM-uuuu HH:mm:ss.SSSSS'
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function rtrn = bcd2char(cBCD)
        % Taken from the Nortek System Integrator 
        % Manual "Example Program" Chapter.
        cBCD = min(cBCD, 153);
        rtrn = bitand(cBCD, 15);
        rtrn = rtrn + 10 * bitshift(cBCD,-4);
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = concat_struct(struct1, struct2)
        out = struct();
        fn1 = fieldnames(struct1);
        fn2 = fieldnames(struct2);
        for i = 1:numel(fn1)
            fni = string(fn1(i));
            out.(fni) = struct1.(fni);
        end
        for i = 1:numel(fn2)
            fni = string(fn2(i));
            out.(fni) = struct2.(fni);
        end
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = fillgaps(a, extrapFlg)
    % Linearly fill NaN value in an array.
    % 
    % Parameters
    % ----------
    % a : array
    %   The array to be filled.    
    % extrapFlg : bool 
    %   Whether to extrapolate if NaNs are found at the ends of the
    %   array.    
    % 
    % Notes
    % -----
    % This function interpolates assuming spacing/timestep between
    % successive points is constant. If the spacing is not constant, use
    % interpgaps.
    % 
        gd = find(~isnan(a));
    
        % Extrapolate if requested
        if extrapFlg && ~isempty(gd)
            if gd(1) ~= 1 
                a(1:gd(1)) = a(gd(1));
            end
            if gd(end) ~= length(a)
                a(gd(end):end) = a(gd(end));
            end
        end
        % Main loop
        if length(gd) > 1
            inds = find(bitand((1 < diff(gd)),(diff(gd)<=inf)));
            for i2 = 1:length(inds)
                ii = gd(inds(i2))+1:gd(inds(i2)+1)-1;
                a(ii) = diff([a(gd(inds(i2))),a(gd(inds(i2)+1))]) *...
                    (1:length(ii)) / (length(ii)) + a(gd(inds(i2)));
            end
        end
        out = a;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function out = interpgaps(a, t, extrapFlg)
    % """
    % Fill gaps (NaN values) in ``a`` by linear interpolation along
    % dimension ``dim`` with the point spacing specified in ``t``.
    % 
    % Parameters
    % ----------
    % a : array
    %   The array containing NaN values to be filled.
    % t : array (len(t) == a.shape[dim])
    %   Independent variable of the points in ``a``, e.g. timestep 
    % extrapFlg : bool 
    %   Whether to extrapolate if NaNs are found at the ends of the
    %   array.
    % 
    % See Also
    % --------
    % dolfyn.io.fillgaps : Linearly interpolates in array-index space.
        
        gd = find(~isnan(a));
    
        % Extrapolate if requested
        if extrapFlg && ~isempty(gd)
            if gd(1) ~= 1 
                a(1:gd(1)) = a(gd(1));
            end
            if gd(end) ~= length(a)
                a(gd(end):end) = a(gd(end));
            end
        end

        % Main loop
        if length(gd) > 1
            inds = find(bitand((1 < diff(gd)),(diff(gd)<=inf)));
            for i2 = 1:length(inds)
                ii = gd(inds(i2))+1:gd(inds(i2)+1)-1;
                ti = (t(ii) - t(gd(inds(i2)))) / ...
                    diff([t(gd(inds(i2))),t(gd(inds(i2)+1))]);
                a(ii) = diff([a(gd(inds(i2))),a(gd(inds(i2)+1))]) *...
                    ti + a(gd(inds(i2)));
            end
        end
        out = a;
    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    function handle_nan()
    % 
    % Finds nan's that cause issues in running the rotation algorithms
    % and deletes them.
    % 
        not_nums = false(size(data.coords.time));
        l = length(data.coords.time);

        if any(isnan(data.coords.time))
            not_nums = not_nums | isnan(data.coords.time);
        end

        var = {'accel', 'angrt', 'mag'};
        fn = fieldnames(data.data_vars);
        for k=1:numel(fn)
            if any(strcmp(fn{k}, var))
                shp = size(data.data_vars.(fn{k}));
                if shp(1) == l
                    if length(shp) == 2
                        if any(isnan(data.data_vars.(fn{k})))
                            temp = (isnan(data.data_vars.(fn{k})))';
                            not_nums = not_nums | temp;                                
                        end
                    elseif length(shp) == 3
                        if any(isnan(data.data_vars.(fn{k})(:,:,end)))
                            temp = isnan(data.data_vars.(fn{k}));
                            temp = squeeze(temp(:,:,end));
                            not_nums = not_nums | temp;
                        end
                    end
                end
            end
        end

        if sum(not_nums,'all') > 0
            data.coords.time = data.coords.time(~not_nums);
            for k=1:numel(fn)
                if length(data.data_vars.(fn{k})) == l
                    dims = size(data.data_vars.(fn{k}));
                    dim = dims(end);
                    dim2 = dims(end-1);
                    if dim < 3
                        data.data_vars.(fn{k}) = ...
                            data.data_vars.(fn{k})(~not_nums);
                    else
                        size_ = size(data.data_vars.(fn{k}));
                        size_(1) = sum(~not_nums,'all');
                        fill = zeros(size_);
                        otherdims = repmat({':'},1,length(dims)-2);
                        for qq = 1:dim                            
                            for jj = 1:dim2
                                temp = data.data_vars.(fn{k})(otherdims{:},jj,qq);
                                fill(otherdims{:},jj,qq) = temp(~not_nums);
                            end
                        end
                        data.data_vars.(fn{k}) = fill;
                    end
                end
            end
        end

    end

    % <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>    

end