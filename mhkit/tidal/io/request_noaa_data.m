function data=request_noaa_data(station, parameter, start_date, end_date, varargin)

%     Loads NOAA current data directly from https://tidesandcurrents.noaa.gov/api/ using a 
%     GET request into a structure. NOAA sets max of 31 days between start and end date.
%     See https://co-ops.nos.noaa.gov/api/ for options. All times are reported as GMT and metric
%     units are returned for data.
% 
%     Parameters
%     ----------
%     station : str
%         NOAA current station number (e.g. 's08010')
%     parameter : str
%         NOAA paramter (e.g. 'currents')
%     start_date : str
%         Start date in the format yyyyMMdd
%     end_date : str
%         End date in the format yyyyMMdd 
%     proxy : structure (optional)
%          To request data from behind a firewall, define a structure of proxy settings, 
%          for example proxy.http = "localhost:8080"
%     write_json : str (optional)
%         Name of json file to write data
%         
%     Returns
%     -------
%     data : structure 
%
%
%         data.id: station ID
%
%         data.name: station name
%
%         data.lat: station Latitude
%
%         data.lon: station Longitude
%
%         data.vars: this will vary depending on parameter input. 
%
%         data.time: epoch time [s]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

if nargin == 5
    if isstruct(varargin{1})
        proxy=py.dict(varargin{1});
        datap=py.mhkit.tidal.io.request_noaa_data(station,parameter,start_date,end_date,pyargs('proxy',proxy));
    elseif isstring(varargin{1})
        datap=py.mhkit.tidal.io.request_noaa_data(station,parameter,start_date,end_date,pyargs('write_json',varargin{1}));
    else
        ME = MException('MATLAB:request_noaa_data','variable argument of wrong data type');
        throw(ME);
    end
elseif nargin == 6
    if isstruct(varargin{1}) & isstring(varargin{2})
        proxy=py.dict(varargin{1});
        datap=py.mhkit.tidal.io.request_noaa_data(station,parameter,start_date,end_date,pyargs('proxy',proxy,'write_json',varargin{2}));
    elseif isstruct(varargin{2}) & isstring(varargin{1})
        proxy=py.dict(varargin{2});
        datap=py.mhkit.tidal.io.request_noaa_data(station,parameter,start_date,end_date,pyargs('proxy',proxy,'write_json',varargin{1}));
    else
        ME = MException('MATLAB:request_noaa_data','wrong vaiable type passed');
        throw(ME);
    end
elseif nargin == 4
    datap=py.mhkit.tidal.io.request_noaa_data(station,parameter,start_date,end_date);
else
    ME = MException('MATLAB:request_noaa_data','incorrect number of variables passed');
        throw(ME);
end
datac=cell(datap);
data=struct(datac{2});
data_df=datac{1};
disp(data_df)

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


