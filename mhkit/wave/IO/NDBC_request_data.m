function ndbc_data=NDBC_request_data(parameter,filenames,options)

%%%%%%%%%%%%%%%%%%%%
%     Requests data by filenames and returns a structure of structures 
%     for each filename passed. If filenames for a single buoy are passed 
%     then the yearly structures in the returned structure (ndbc_data) are 
%     indexed by year (e.g. ndbc_data.year_2004). If multiple buoy ids are 
%     passed then the returned dictionary is indexed by buoy id and year 
%     (e.g. ndbc_data['46022']['2014']). 
%     
%     
% Parameters
% ------------
%     parameter : string
%         'swden'	:	'Raw Spectral Wave Current Year Historical Data'
%         'stdmet':   'Standard Meteorological Current Year Historical Data'
%
%     filenames : array of strings
%         Data filenames on https://www.ndbc.noaa.gov/data/historical/{parameter}/  
%         
%
%     proxy : structure or py.None (optional)
%         To request data from behind a firewall, define a dictionary of proxy settings, 
%         for example proxy.http = "http:wwwproxy.yourProxy:80/"
%         to call: NDBC_request_data(parameter,filenames,"proxy",proxy)
%     
% Returns
% ---------
%     ndbc_data: Structure 
%         Structure of structures broken down by years of data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    parameter
    filenames
    options.proxy = py.None;
end

py.importlib.import_module('mhkit');

if (isa(options.proxy,'py.NoneType')~=1)
    options.proxy=py.dict(options.proxy);
end

li=py.list();

 ind = 1:1:length(filenames);
 for k=1:length(filenames)
     li=py.mhkit_python_utils.pandas_dataframe.lis(li,filenames(k));
 end

df = py.pandas.DataFrame(pyargs("data",li));

datapd = py.mhkit.wave.io.ndbc.request_data(parameter,df);
key = cell(py.list(datapd));

if (isa(datapd{string(py.str(key{1}))}.values,'py.dict_values')==1)
   for i=1:length(key)
    a = char(string(py.str(key{i})));
    dict2 = datapd{string(py.str(key{i}))};
    key2 = cell(py.list(dict2));
    for j=1:length(key2)
       df = double(dict2{string(py.str(key2{j}))}.values);
       cols = cell(py.list(py.numpy.nditer(dict2{string(py.str(key2{j}))}.columns,pyargs("flags",{"refs_ok"})))); 
       for k = 1:length(cols)
        b = char(string(py.str(key2{j})));
        c = char(string(py.str(cols{k})));
        d = df(:,k);
        eval(['ndbc_data.ID_' a '.year_' b '.' c '= d ;']);
       end
    end
   end
 
else
    for i=1:length(key)
        df = double(datapd{string(py.str(key{i}))}.values);
        cols = cell(py.list(py.numpy.nditer(datapd{string(py.str(key{i}))}.columns,pyargs("flags",{"refs_ok"}))));
        for k = 1:length(cols)
            a = char(string(py.str(key{i})));
            b = char(string(py.str(cols{k})));
            c = df(:,k);
            eval(['ndbc_data.year_' a '.' b '= c ;']);
        end
    end
end
