function ndbc_data=NDBC_request_data(parameter,filenames,options)

%%%%%%%%%%%%%%%%%%%%
%     Requests data by filenames and returns a dictionary of DataFrames 
%     for each filename passed. If filenames for a sigle buoy are passed 
%     then the yearly DataFrames in the returned dictionary (ndbc_data) are 
%     indexed by year (e.g. ndbc_data['2014']). If multiple buoy ids are 
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
%         Buoy Number.  5-character alpha-numeric station identifier 
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
%         Structure with station ID, years, and NDBC file names.
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
%key = cell(py.list(py.numpy.nditer(keys(datapd),pyargs("flags",{"refs_ok"}))));
disp(string(py.str(key{1})))
ndbc_data = 2;