function ds = create_dataset(data)
    % Creates final dataset from structure created from binary readers.
    % It is meant to try and approximate an xarray dataset
    % Direction 'dir' coordinates get reset in `set_coords`
    ds = struct;   
    inst = {'X' 'Y' 'Z'};
    earth = {'E' 'N' 'U'};
    beam = 1:length(size(data.data_vars.vel));
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
                ds.(fn{k}).dims = {'x', 'x_star'};
                ds.(fn{k}).coords.x = beam;
                ds.(fn{k}).coords.x_star = beam;
            else 
            % earth2inst orientation matrix  
                ds.(fn{k}).data = data.data_vars.(fn{k});
                ds.(fn{k}).dims = {strcat('time',tg), 'earth', 'inst' };
                ds.(fn{k}).coords.(strcat('time',tg)) = ...
                    data.coords.(strcat('time',tg));
                ds.(fn{k}).coords.earth = earth;
                ds.(fn{k}).coords.inst = inst;                
            end
            
        % quaternion units never change
        elseif (contains(fn{k},'quaternion'))
            
            ds.(fn{k}).data = data.data_vars.(fn{k});
            ds.(fn{k}).dims = {strcat('time',tg), 'q'};
            ds.(fn{k}).coords.(strcat('time',tg))= ...
                    data.coords.(strcat('time',tg));
            ds.(fn{k}).coords.q = {'w', 'x', 'y', 'z'};            
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

