function datast=read_NDBC_file(file_name,varargin)

%%%%%%%%%%%%%%%%%%%%
%     Reads a NDBC wave buoy data file (from https://www.ndbc.noaa.gov) into a
%     structure.
%     
%     Realtime and historical data files can be loaded with this function.  
%     
%     Note: With realtime data, missing data is denoted by "MM".  With historical 
%     data, missing data is denoted using a variable number of 
%     # 9's, depending on the data type (for example: 9999.0 999.0 99.0).
%     'N/A' is automatically converted to missing data.
%     
%     Data values are converted to float/int when possible. Column names are 
%     also converted to float/int when possible (this is useful when column 
%     names are frequency).
%     
% Parameters
% ------------
%     file_name : string
%         Name of NDBC wave buoy data file
%
%     missing_value : vector of values (optional)
%         vector of values that denote missing data    
%     
% Returns
% ---------
%     data: Structure 
%
%
%         data.Data: named according to header row 
%
%         data.time: given in datetime 
%
%         data.units: the units for each data entry
%         
%         OR if a spectra data NDBC file
%
%         data.spectra: spectra data
%
%         data.time: given in datetime 
%
%         data.frequency: spectral frequency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

if nargin == 2
    missing=py.list(varargin{1});
    datatp=py.mhkit.wave.io.ndbc.read_file(file_name,missing);
elseif nargin > 2
    ME = MException('MATLAB:read_NDBC_file','too many arguments passed');
    throw(ME);
else
    datatp=py.mhkit.wave.io.ndbc.read_file(file_name);
end

datac=cell(datatp);
datapd=datac{1};

datamat=datac{2};
matstr=struct(datamat);

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
temp = [];
if ~isempty(fieldnames(matstr))
    for i=1:si(2)
        test=string(py.str(vv{i}));
        datast.(test)=vals(:,i);
        unit=string(matstr.(test));
        datast.units.(test)=unit;
    end
else
    datast.spectrum = vals';
    for i=1:si(2)
        temp = [temp, double(py.array.array('d',py.numpy.nditer(vv{i})))];
    end
    datast.frequency = temp';
end

 for i=1:siti(2)
    datast.time(i)=datetime(string(py.str(ti{i})),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS')';
 end

