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

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

datap = py.mhkit.tidal.io.read_noaa_json(filename);

datac=cell(datap);
data=struct(datac{2});
data_df=datac{1};

fields=fieldnames(data);
for idx = 1:length(fields)
    data.(fields{idx}) = string(data.(fields{idx}));
end

xx=cell(data_df.axes);
v=xx{2};


vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(data_df.values,pyargs("flags",{"refs_ok"}))));
sha=cell(data_df.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);

ti=cell(py.list(py.numpy.nditer(data_df.index,pyargs("flags",{"refs_ok"}))));
siti=size(ti);
si=size(vals);


 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");
    
    data.(newname(1))=vals(:,i);
    
 end
 for i=1:siti(2)
    data.time{i}=posixtime(datetime(string(py.str(ti{i})),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS'));
 end

data.time=cell2mat(data.time);