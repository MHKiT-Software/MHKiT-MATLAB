function datast=request_usgs_data(station, parameter, start_date, end_date,varargin)

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
%
%     proxy : dict or None (optional)
%          To request data from behind a firewall, define a dictionary of proxy settings, 
%          for example {"http": 'localhost:8080'}
%
%     write_json : str or None (optional)
%         Name of json file to write data
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

%%%%% NOTE: TODO: still need to add creation of dictionary from structure
%%%%% for proxy argument 

py.importlib.import_module('mhkit');

if nargin > 4
    if nargin > 7
        ME = MException('MATLAB:request_usgs_data','Incorrect number of input arguments, too many agruments, requires 7 at most, %d arguments passed',nargin);
        throw(ME);
    
    elseif nargin == 5
        
            if any([strcmp(varargin{1},'Daily'), strcmp(varargin{1},'Instantaneous')]) 
                datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{1}));
            elseif contains(varargin{1},'.json')
                datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('write_json',varargin{1}));
            
        elseif (isa(varargin{1},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('proxy',varargin{1}));
        else
           ME = MException('MATLAB:request_usgs_data','One or more optional argument is of the wrong type');
           throw(ME); 
        end
        
        
    elseif nargin == 6
        if any([strcmp(varargin{1},'Daily'), strcmp(varargin{1},'Instantaneous')]) & contains(varargin{2},'.json')
        	datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{1}, 'write_json',varargin{2}));
        elseif any([strcmp(varargin{2},'Daily'), strcmp(varargin{2},'Instantaneous')]) & contains(varargin{1},'.json')
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{2}, 'write_json',varargin{1}));
        elseif any([strcmp(varargin{1},'Daily'), strcmp(varargin{1},'Instantaneous')]) & (isa(varargin{2},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{1}, 'proxy',varargin{2}));
        elseif any([strcmp(varargin{2},'Daily'), strcmp(varargin{2},'Instantaneous')]) & (isa(varargin{1},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{2}, 'proxy',varargin{1}));
        elseif (isa(varargin{1},'py.dict')==1) & contains(varargin{2},'.json')
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('write_json',varargin{2}, 'proxy',varargin{1}));
        elseif (isa(varargin{2},'py.dict')==1) & contains(varargin{1},'.json')
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('write_json',varargin{1}, 'proxy',varargin{2}));
        else 
            ME = MException('MATLAB:request_usgs_data','One or more optional argument is of the wrong type');
            throw(ME);
        end
    elseif nargin == 7 
        if any([strcmp(varargin{1},'Daily'), strcmp(varargin{1},'Instantaneous')]) & contains(varargin{2},'.json') & (isa(varargin{3},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{1},...
                'write_json',varargin{2}, 'proxy',varargin{3}));
        elseif any([strcmp(varargin{1},'Daily'), strcmp(varargin{1},'Instantaneous')]) & contains(varargin{3},'.json') & (isa(varargin{2},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{1},...
                'write_json',varargin{3}, 'proxy',varargin{2}));
        elseif any([strcmp(varargin{2},'Daily'), strcmp(varargin{2},'Instantaneous')]) & contains(varargin{3},'.json') & (isa(varargin{1},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{2},...
                'write_json',varargin{3}, 'proxy',varargin{1}));
        elseif any([strcmp(varargin{2},'Daily'), strcmp(varargin{2},'Instantaneous')]) & contains(varargin{1},'.json') & (isa(varargin{3},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{2},...
                'write_json',varargin{1}, 'proxy',varargin{3}));
        elseif any([strcmp(varargin{3},'Daily'), strcmp(varargin{3},'Instantaneous')]) & contains(varargin{1},'.json') & (isa(varargin{2},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{3},...
                'write_json',varargin{1}, 'proxy',varargin{2}));
        elseif any([strcmp(varargin{3},'Daily'), strcmp(varargin{3},'Instantaneous')]) & contains(varargin{2},'.json') & (isa(varargin{1},'py.dict')==1)
            datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date,pyargs('data_type',varargin{3},...
                'write_json',varargin{2}, 'proxy',varargin{1}));
        else
            ME = MException('MATLAB:request_usgs_data','One or more optional argument is of the wrong type');
            throw(ME);
        end
    end
    
else
    datapd=py.mhkit.river.io.request_usgs_data(station, parameter, start_date, end_date);
end


xx=cell(datapd.axes);
v=xx{2};


vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(datapd.values)));
sha=cell(datapd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);
ti=cell(py.list(py.numpy.nditer(datapd.index)));
siti=size(ti);
si=size(vals);
 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");
    
    datast.(newname(1))=vals(:,i);
    datast.units.(newname(1))=newname(2);
 end
 for i=1:siti(2)
    datast.time{i}=posixtime(datetime(string(py.str(ti{i})),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS'));
 end

datast.time=cell2mat(datast.time);

