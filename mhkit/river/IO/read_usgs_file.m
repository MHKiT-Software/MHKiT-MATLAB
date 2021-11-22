function datast=read_usgs_file(file_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reads a USGS JSON data file (from https://waterdata.usgs.gov/nwis) 
%     into a structure 
% 
% Parameters
% ----------
%     file_name : str
%         Name of USGS JSON data file
%         
% Returns
% -------
%     datast : structure
%
%
%         datast.Data: named according to the parameter's variable description
%
%         datast.time: epoch time [s]
%
%         datast.units: units for each parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');

datapd=py.mhkit.river.io.usgs.read_usgs_file(file_name);


xx=cell(datapd.axes);
v=xx{2};


vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(datapd.values,pyargs("flags",{"refs_ok"}))));
sha=cell(datapd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);
ti=cell(py.list(py.numpy.nditer(datapd.index,pyargs("flags",{"refs_ok"}))));
siti=size(ti);
si=size(vals);
 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");
    
    datast.(newname(1))=vals(:,i);
    datast.units.(newname(1))=newname(2);
 end
 for i=1:siti(2)
    datast.time{i}=posixtime(datetime(string(py.str(ti{i})),'InputFormat','yyyy-MM-dd HH:mm:ssXXX','TimeZone','UTC'));
 end
 
 datast.time=cell2mat(datast.time);



