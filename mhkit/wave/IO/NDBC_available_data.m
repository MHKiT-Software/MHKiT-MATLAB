function available_data=NDBC_available_data(parameter,options)

%%%%%%%%%%%%%%%%%%%%
%     For a given parameter this will return a structure of years, 
%     station IDs and file names that contain that parameter data. 
%     
%     
% Parameters
% ------------
%     parameter : string
%         'swden'	:	'Raw Spectral Wave Current Year Historical Data'
%         'stdmet':   'Standard Meteorological Current Year Historical Data'
%
%     buoy_number : string or py.None (optional)
%         Buoy Number.  5-character alpha-numeric station identifier 
%         to call: NDBC_available_data(parameter,"buoy_number",buoy_number)
%
%     proxy : structure or py.None (optional)
%         To request data from behind a firewall, define a dictionary of proxy settings, 
%         for example proxy.http = "http:wwwproxy.yourProxy:80/"
%         to call: NDBC_available_data(parameter,"proxy",proxy)
%     
% Returns
% ---------
%     avaialble_data: Structure 
%         Structure with station ID, years, and NDBC file names.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    parameter
    options.buoy_number = py.None;
    options.proxy = py.None;
end

py.importlib.import_module('mhkit');

if (isa(options.proxy,'py.NoneType')~=1)
    options.proxy=py.dict(options.proxy);
end

datapd=py.mhkit.wave.io.ndbc.available_data(parameter,pyargs('buoy_number',options.buoy_number,...
      'proxy',options.proxy));
datali = cell(py.list(py.numpy.nditer(datapd.values,pyargs("flags",{"refs_ok"}))));
sha=cell(datapd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});


siti=size(datali);
for i=1:siti(2)
    st{i} = string(py.str(datali{i}));

end
st2 = reshape(st,x,y);

station_id = [st2{:,1}];
year = [st2{:,2}];
file = [st2{:,3}];

varnames = {'Station_id','year','file'};
available_data = table(station_id',year',file','VariableNames',varnames);