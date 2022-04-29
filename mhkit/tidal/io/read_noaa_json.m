function data=read_noaa_json(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns site structure from a json saved from the 
%     request_noaa_data
%
% Parameters
% ----------
%     filename: string
%         filename with path of json file to load
% Returns
% -------
%     data: Structure
%
%
%         data.Data: Timeseries Site data, will be named based on parameters in JSON file
%
%         data.time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');
% 
% datap = py.mhkit.tidal.io.noaa.read_noaa_json(filename);
% 
% datac=cell(datap);
% data=struct(datac{2});
% data_df=datac{1};
% 
% fields=fieldnames(data);
% 
% for idx = 1:length(fields)
%     data.(fields{idx}) = string(data.(fields{idx}));
% end
% 
% xx=cell(data_df.axes);
% v=xx{2};
% 
% 
% vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));
% 
% vals=double(py.array.array('d',py.numpy.nditer(data_df.values,pyargs("flags",{"refs_ok"}))));
% sha=cell(data_df.values.shape);
% x=int64(sha{1,1});
% y=int64(sha{1,2});
% 
% vals=reshape(vals,[x,y]);
% si=size(vals);
% 
% for i=1:si(2)
%     test=string(py.str(vv{i}));
%     newname=split(test,",");
%     data.(newname(1))=vals(:,i);
% end
%  
% times = double(    ...
%      py.mhkit_python_utils.pandas_dataframe.datetime_index_to_ordinal(data_df));
% 
% data.time = posixtime(datetime(times, 'ConvertFrom', 'datenum'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);  % Open the file
raw = fread(fid,inf);   % Read the contents
str = char(raw');       % Transformation
fclose(fid);            % Close the file
init_data = jsondecode(str);

data = struct();

% Move metadata from its own field to the main structure
fields = fieldnames(init_data.metadata);
for idx = 1:length(fields)
    data.(fields{idx}) = convertCharsToStrings(...
        init_data.metadata.(fields{idx}));
end
init_data = rmfield(init_data,'metadata');

% Move the data fields over and create time 
fields = fieldnames(init_data);
shape = size(fieldnames(init_data.(fields{1})));
for idx = 1:length(fields)
    data.(fields{idx}) = zeros(shape);
end

% loop through the time (which gets read in as the field names for each
% data column) and fill in data

time = fieldnames(init_data.(fields{1}));
data.time = cellfun(@(x) x(2:end), time, 'UniformOutput', false);
data.time = cellfun(@str2double,data.time);
data.time = transpose(data.time/1000.);
% Time is in epochtime since 1970-01-01. Divide by 1000 to remove
% milliseconds. convert with 
% datetime(time,'ConvertFrom','epochtime','Epoch','1970-01-01')

for idx = 1:length(fields)
    test = init_data.(fields{idx}).(time{1});    
    data.(fields{idx}) = cell2mat(struct2cell(init_data.(fields{idx})));
    if ischar(test) || isstring(test)        
        data.(fields{idx}) = arrayfun(@str2double,data.(fields{idx}));
    end
end

end