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