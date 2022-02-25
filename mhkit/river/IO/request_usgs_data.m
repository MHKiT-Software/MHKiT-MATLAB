function datast=request_usgs_data(station, parameter, start_date, end_date,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Loads USGS data directly from https://waterdata.usgs.gov/nwis into a structure using a 
%     GET request 
%     
% 
% Parameters
% ----------
%     station : str
%         USGS station number (e.g. '08313000')
%
%     parameter : str
%         USGS paramter ID (e.g. '00060' for Discharge, cubic feet per second)
%
%     start_date : str
%         Start date in the format 'YYYY-MM-DD' (e.g. '2018-01-01')
%
%     end_date : str
%         End date in the format 'YYYY-MM-DD' (e.g. '2018-12-31')
%
%     data_type : str (optional)
%         Data type, options include 'Daily' (return the mean daily value) and 
%         'Instantaneous'.
%         to call: request_usgs_data(station,parameter,start_date.end_date,"data_type",data_type)
%
%     proxy : structure or None (optional)
%         To request data from behind a firewall, define a dictionary of proxy settings, 
%         for example proxy.http = "localhost:8080"
%         to call: request_usgs_data(station,parameter,start_date.end_date,"proxy",proxy)
%
%     write_json : str or None (optional)
%         Name of json file to write data
%         to call: request_usgs_data(station,parameter,start_date.end_date,"write_json",write_json)
%         
% Returns
% -------
%     datast : structure 
%
%
%         datast.Data: named according to the parameter's variable description
%
%         datast.time: datetime
%
%         datast.units: units for each parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    station 
    parameter
    start_date
    end_date
    options.data_type = 'Daily';
    options.proxy = py.None;
    options.write_json = py.None; 
end

py.importlib.import_module('mhkit');

if (isa(options.proxy,'py.NoneType')~=1)
    options.proxy=py.dict(options.proxy);
end
datapd=py.mhkit.river.io.usgs.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',options.data_type,...
      'write_json',options.write_json, 'proxy',options.proxy));

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

