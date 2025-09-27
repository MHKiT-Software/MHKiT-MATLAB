function ds = create_dataset(data)
    % Creates final dataset from structure created from binary readers.
    % It is meant to try and approximate an xarray dataset
    % Direction 'dir' coordinates get reset in `set_coords`
    ds = struct;
    inst = {'X' 'Y' 'Z'};
    earth = {'E' 'N' 'U'};
    beam = 1:min(size(squeeze(data.data_vars.vel)));
    tag = {'_b5', '_echo', '_bt', '_gps', '_ast'};
    ds.coords.inst = inst;
    ds.coords.earth = earth;

    fn = fieldnames(data.data_vars);
    for k=1:numel(fn)
        tmp = split(fn{k},"_");
        chk_tag = strcat("_",tmp(end));
        if any(strcmp(tag,chk_tag))
            tg = chk_tag;
        else
            tg = '';
        end
        if (contains(fn{k},'mat'))
            % beam2inst & inst2head orientation matrices
            if (contains(fn{k},'inst'))
                ds.(fn{k}).data = data.data_vars.(fn{k});
                ds.(fn{k}).dims = {'x1', 'x2'};
                ds.(fn{k}).coords.x1 = beam;
                ds.(fn{k}).coords.x2 = beam;
                % Add units for beam2inst_orientmat
                ds.(fn{k}).units = '1';
            else
            % earth2inst orientation matrix
                ds.(fn{k}).data = data.data_vars.(fn{k});
                ds.(fn{k}).dims = {strcat('time',tg), 'inst', 'earth'};
                ds.(fn{k}).coords.(strcat('time',tg)) = ...
                    data.coords.(strcat('time',tg));
                ds.(fn{k}).coords.inst = inst;
                ds.(fn{k}).coords.earth = earth;
                % Add units for orientmat
                ds.(fn{k}).units = '1';
            end

        % quaternion units never change
        elseif (contains(fn{k},'quaternion'))

            ds.(fn{k}).data = data.data_vars.(fn{k});
            ds.(fn{k}).dims = {strcat('time',tg), 'q'};
            ds.(fn{k}).coords.(strcat('time',tg))= ...
                    data.coords.(strcat('time',tg));
            ds.(fn{k}).coords.q = {'w', 'x', 'y', 'z'};

        % Handle altraw variables specially to add time_altraw coordinates
        elseif contains(fn{k}, 'altraw') && ~contains(fn{k}, 'samp_altraw')
            ds.(fn{k}).data = data.data_vars.(fn{k});
            if isfield(data.units,fn{k})
                ds.(fn{k}).units = data.units.(fn{k});
            else
                % Special handling for altraw pressure
                if contains(fn{k}, 'pressure_altraw')
                    ds.(fn{k}).units = 'dbar';
                end
            end

            % Determine tag for altraw coordinates
            if contains(fn{k}, '_avg')
                altraw_tag = '_avg';
            else
                altraw_tag = '';
            end

            % Add time_altraw coordinate
            time_altraw_key = strcat('time_altraw', altraw_tag);
            if isfield(data.coords, time_altraw_key)
                shp = size(data.data_vars.(fn{k}));
                if length(shp) == 2 && shp(2) == 1
                    % 1D altraw variable
                    ds.(fn{k}).dims = {time_altraw_key};
                    ds.(fn{k}).coords.(time_altraw_key) = data.coords.(time_altraw_key);
                end
            end

        % Handle altimeter data fields specially
        elseif contains(fn{k}, '_alt') && ~contains(fn{k}, 'altraw')
            ds.(fn{k}).data = data.data_vars.(fn{k});
            if isfield(data.units,fn{k})
                ds.(fn{k}).units = data.units.(fn{k});
            else
                % Set default units for known altimeter fields
                if contains(fn{k}, 'pressure_alt')
                    ds.(fn{k}).units = 'dbar';
                elseif contains(fn{k}, 'dist_alt')
                    ds.(fn{k}).units = 'm';
                elseif contains(fn{k}, 'quality_alt')
                    ds.(fn{k}).units = 'dB';
                elseif contains(fn{k}, 'offset_time_alt')
                    ds.(fn{k}).units = 's';
                else
                    ds.(fn{k}).units = '';
                end
            end

            % Determine tag and time coordinate for altimeter data
            if contains(fn{k}, '_avg')
                alt_tag = '_avg';
                time_coord = 'time_avg';
            else
                alt_tag = '';
                time_coord = 'time';
            end

            % Add time coordinate for 1D altimeter variables
            shp = size(data.data_vars.(fn{k}));
            if length(shp) == 2 && shp(2) == 1
                % 1D altimeter variable - use regular time coordinate
                if isfield(data.coords, time_coord)
                    ds.(fn{k}).dims = {time_coord};
                    ds.(fn{k}).coords.(time_coord) = data.coords.(time_coord);
                end
            end
        else
            ds.(fn{k}).data = data.data_vars.(fn{k});
            if isfield(data.units,fn{k})
                ds.(fn{k}).units = data.units.(fn{k});
            else
                % make sure ones with tags get units
                if ~isempty(tg)
                    ds.(fn{k}).units =...
                        data.units.(strjoin([tmp(1:end-1)],'_'));
                end
            end

            shp = size(data.data_vars.(fn{k}));
            vshp = size(data.data_vars.vel);
            l = length(shp);

            if l == 2 % 1D variables
                if any(strcmp(tag,chk_tag))
                    tg = chk_tag;
                else
                    tg = '';
                end
                % I'm not sure what this part does. Need to come back
                % to it if its used for the other data types
                ds.(fn{k}).dims = {strcat("time",tg)};
                ds.(fn{k}).coords.(strcat('time',tg)) = ...
                    data.coords.(strcat('time',tg));

            elseif l == 3 % 2D variables
                sub_tag = tag(1:2);
                if strcmp('echo',fn{k})
                    ds.(fn{k}).dims = {'time_echo', 'range_echo'};
                    ds.(fn{k}).coords.time_echo = ...
                        data.coords.time_echo;
                    ds.(fn{k}).coords.range_echo = ...
                        data.coords.range_echo;

                % 3- & 4-beam instrument vector data, bottom tracking
                elseif shp(end)==vshp(end) && ~any(strcmp(sub_tag,chk_tag))
                    % b/c rdi time
                    if (contains(fn{k},'bt')) && ...
                        isfield(data.coords,'time_bt')
                        tg = '_bt';
                    else
                        tg = '';
                    end
                    ds.(fn{k}).dims = {strcat("time",tg), 'dir'};
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));
                    ds.(fn{k}).coords.dir = beam;
                % 4-beam instrument IMU data
                elseif shp(end) == vshp(end)-1
                    ds.(fn{k}).dims = {strcat("time",tg), 'dirIMU'};
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));
                    ds.(fn{k}).coords.dirIMU = 1:3;

                elseif any(strcmp(sub_tag,chk_tag))
                    ds.(fn{k}).dims = {strcat("time",tg),...
                        strcat("range",tg)};
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));
                    ds.(fn{k}).coords.(strcat("range",tg)) = ...
                        data.coords.(strcat('range',tg));
                else
                    warning(strcat('Variable not included in dataset: ',...
                        fn{k}));
                end

            elseif l == 4
               if ~contains(fn{k},join(tag))
                   if (contains(fn{k},'vel'))
                       dim0 = 'dir';
                   else
                       dim0 = 'beam';
                   end
                   ds.(fn{k}).dims = {'time','range', dim0};
                   ds.(fn{k}).coords.time = data.coords.time;
                   ds.(fn{k}).coords.range = data.coords.range;
                   ds.(fn{k}).coords.(dim0) = beam;

               elseif (contains(fn{k},'b5'))
                   ds.(fn{k}).dims = {'time_b5', 'range_b5' };
                   ds.(fn{k}).coords.time_b5 = data.coords.time_b5;
                   ds.(fn{k}).coords.range_b5 = data.coords.range_b5;
               else
                   warning(strcat('Variable not included in dataset: ',...
                        fn{k}));
               end
            end
        end
        add_to_coords(fn{k});
    end

    % coordinate units
    ds.coord_sys = data.attrs.coord_sys;

    ds.attrs = data.attrs;
    ds.time = data.coords.time;
    if isfield(ds,'range')
        ds.range = data.coords.range;
    end
    if isfield(data.coords,'n_altraw')
        ds.n_altraw = data.coords.n_altraw;
        ds.coords.n_altraw = data.coords.n_altraw;
    end

    % Add sampling structures for altraw data if they exist
    % nsamp_altraw structure
    if isfield(data.data_vars, 'nsamp_altraw')
        ds.nsamp_altraw.data = data.data_vars.nsamp_altraw;
        ds.nsamp_altraw.dims = {'time_altraw'};
        if isfield(data.coords, 'time_altraw')
            ds.nsamp_altraw.coords.time_altraw = data.coords.time_altraw;
        end
        if isfield(data.units, 'nsamp_altraw')
            ds.nsamp_altraw.units = data.units.nsamp_altraw;
        else
            ds.nsamp_altraw.units = '';
        end
    end

    % dsamp_altraw structure
    if isfield(data.data_vars, 'dsamp_altraw')
        ds.dsamp_altraw.data = data.data_vars.dsamp_altraw;
        ds.dsamp_altraw.dims = {'time_altraw'};
        if isfield(data.coords, 'time_altraw')
            ds.dsamp_altraw.coords.time_altraw = data.coords.time_altraw;
        end
        if isfield(data.units, 'dsamp_altraw')
            ds.dsamp_altraw.units = data.units.dsamp_altraw;
        else
            ds.dsamp_altraw.units = 'm';
        end
    end

    % samp_altraw structure
    if isfield(data.data_vars, 'samp_altraw')
        ds.samp_altraw.data = data.data_vars.samp_altraw;
        ds.samp_altraw.dims = {'n_altraw', 'time_altraw'};
        if isfield(data.coords, 'n_altraw')
            ds.samp_altraw.coords.n_altraw = data.coords.n_altraw;
        end
        if isfield(data.coords, 'time_altraw')
            ds.samp_altraw.coords.time_altraw = data.coords.time_altraw;
        end
        if isfield(data.units, 'samp_altraw')
            ds.samp_altraw.units = data.units.samp_altraw;
        else
            ds.samp_altraw.units = '';
        end
    end

    % Handle _avg versions as well
    if isfield(data.data_vars, 'nsamp_altraw_avg')
        ds.nsamp_altraw_avg.data = data.data_vars.nsamp_altraw_avg;
        ds.nsamp_altraw_avg.dims = {'time_altraw_avg'};
        if isfield(data.coords, 'time_altraw_avg')
            ds.nsamp_altraw_avg.coords.time_altraw_avg = data.coords.time_altraw_avg;
        end
        ds.nsamp_altraw_avg.units = '';
    end

    if isfield(data.data_vars, 'dsamp_altraw_avg')
        ds.dsamp_altraw_avg.data = data.data_vars.dsamp_altraw_avg;
        ds.dsamp_altraw_avg.dims = {'time_altraw_avg'};
        if isfield(data.coords, 'time_altraw_avg')
            ds.dsamp_altraw_avg.coords.time_altraw_avg = data.coords.time_altraw_avg;
        end
        ds.dsamp_altraw_avg.units = 'm';
    end

    if isfield(data.data_vars, 'samp_altraw_avg')
        ds.samp_altraw_avg.data = data.data_vars.samp_altraw_avg;
        ds.samp_altraw_avg.dims = {'n_altraw', 'time_altraw_avg'};
        if isfield(data.coords, 'n_altraw')
            ds.samp_altraw_avg.coords.n_altraw = data.coords.n_altraw;
        end
        if isfield(data.coords, 'time_altraw_avg')
            ds.samp_altraw_avg.coords.time_altraw_avg = data.coords.time_altraw_avg;
        end
        ds.samp_altraw_avg.units = '';
    end

    function add_to_coords(key)
        fields = fieldnames(ds.(key).coords);
        for qq = 1:numel(fields)
            if ~isfield(ds.coords, fields{qq})
                ds.coords.(fields{qq}) = ds.(key).coords.(fields{qq});
            end
        end
    end

end

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

