function ds = create_dataset(data)
    % Creates final dataset from structure created from binary readers.
    % It is meant to try and approximate an xarray dataset
    % Direction 'dir' coordinates get reset in `set_coords`
    ds = struct;   
    inst = {'X' 'Y' 'Z'};
    earth = {'E' 'N' 'U'};
    beam = 1:length(size(data.data_vars.vel));
    tag = {'_b5', '_echo', '_bt', '_gps', '_ast'};
    
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
                ds.(fn{k}).dims = {'beam' 'x'};
                ds.(fn{k}).coords.beam = beam;
                ds.(fn{k}).coords.x = beam;
            else 
            % earth2inst orientation matrix  
                ds.(fn{k}).data = data.data_vars.(fn{k});
                ds.(fn{k}).dims = {'earth', 'inst', strcat('time',tg)};
                ds.(fn{k}).coords.earth = earth;
                ds.(fn{k}).coords.inst = inst;
                ds.(fn{k}).coords.(strcat('time',tg)) = ...
                    data.coords.(strcat('time',tg));
            end

        % quaternion units never change
        elseif (contains(fn{k},'quaternion'))
            
            ds.(fn{k}).data = data.data_vars.(fn{k});
            ds.(fn{k}).dims = {'q', strcat('time',tg)};
            ds.(fn{k}).coords.q = {'w', 'x', 'y', 'z'};
            ds.(fn{k}).coords.(strcat('time',tg))= ...
                    data.coords.(strcat('time',tg));
        else
            ds.(fn{k}).data = data.data_vars.(fn{k});
            if isfield(data.units,fn{k})
                ds.(fn{k}).units = data.units.(fn{k});
            else
                % make sure ones with tags get units
                if ~isempty(tg) 
                    ds.(fn{k}).units = data.units.(tmp{1});
                else
                    continue
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
                    ds.(fn{k}).dims = {'range_echo' 'time_echo'};
                    ds.(fn{k}).coords.range_echo = ...
                        data.coords.range_echo;
                    ds.(fn{k}).coords.time_echo = ...
                        data.coords.time_echo;

                % 3- & 4-beam instrument vector data, bottom tracking                
                elseif shp(end)==vshp(end) && ~any(strcmp(sub_tag,chk_tag))
                    % b/c rdi time
                    if (contains(fn{k},'bt')) && ...
                        isfield(data.coords,'time_bt')
                        tg = '_bt';
                    else
                        tg = '';
                    end
                    ds.(fn{k}).dims = {'dir', strcat("time",tg)};
                    ds.(fn{k}).coords.dir = beam;
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));
                % 4-beam instrument IMU data
                elseif shp(end) == vshp(end)-1                    
                    ds.(fn{k}).dims = {'dirIMU', strcat("time",tg)};
                    ds.(fn{k}).coords.dirIMU = 1:3;
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));

                elseif any(strcmp(sub_tag,chk_tag))                    
                    ds.(fn{k}).dims = {strcat("range",tg),...
                        strcat("time",tg)};
                    ds.(fn{k}).coords.(strcat("range",tg)) = ...
                        data.coords.(strcat('range',tg));
                    ds.(fn{k}).coords.(strcat("time",tg)) = ...
                        data.coords.(strcat('time',tg));
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
                   ds.(fn{k}).dims = {dim0, 'range', 'time'};
                   ds.(fn{k}).coords.(dim0) = beam;
                   ds.(fn{k}).coords.range = data.coords.range;
                   ds.(fn{k}).coords.time = data.coords.time;

               elseif (contains(fn{k},'b5')) 
                   ds.(fn{k}).dims = {'range_b5', 'time_b5'};
                   ds.(fn{k}).coords.range_b5 = data.coords.range_b5;
                   ds.(fn{k}).coords.time_b5 = data.coords.time_b5;
               else
                   warning(strcat('Variable not included in dataset: ',...
                        fn{k}));
               end
            end
        end  
    end

    % coordinate units
    coords = {};
    fn = fieldnames(ds);
    for k=1:numel(fn)
        coords = [coords ; fieldnames(ds.(fn{k}).coords)];
    end
    coords = unique(coords);
    ds.coords = coords;
    ds.coord_sys = data.attrs.coord_sys;

    ds.attrs = data.attrs;
    ds.time = data.coords.time;
    if isfield(ds,'range')
        ds.range = data.coords.range;
    end

end

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

