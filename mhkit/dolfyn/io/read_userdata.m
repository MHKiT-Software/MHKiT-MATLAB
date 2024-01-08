% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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
            basefile = extractBefore(filename,...
                find(filename == '.', 1, 'last'));
            jsonfile = basefile + ".userdata.json";
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
