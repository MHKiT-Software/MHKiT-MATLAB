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
%         Structure of structures broken down by buoy and years of data.
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
freq = [];
spectra = [];
if (isa(datapd{string(py.str(key{1}))}.values,'py.dict_values')==1)
   for i=1:length(key)
    a = char(string(py.str(key{i})));
    dict2 = datapd{string(py.str(key{i}))};
    key2 = cell(py.list(dict2));
    for j=1:length(key2)
        sha=cell(dict2{string(py.str(key2{j}))}.values.shape);
        x=int64(sha{1,1});
        y=int64(sha{1,2});
       df = reshape(double(py.array.array('d',py.numpy.nditer(dict2{string(py.str(key2{j}))}.values,...
            pyargs("flags",{"refs_ok"})))),[x,y]);
       cols = cell(py.list(py.numpy.nditer(dict2{string(py.str(key2{j}))}.columns,pyargs("flags",{"refs_ok"})))); 
       for k = 1:length(cols)
        
        b = char(string(py.str(key2{j})));
        c = char(string(py.str(cols{k})));
        d = df(:,k);
        if parameter == "swden"
                
                if ~contains(c, ".")
                    if contains(c, "YY")
                        c = char("YYYY");
                    end
                    eval(['ndbc_data.ID_' a '.year_' b '.' c '= d ;']);
                else
                    freq = [freq str2double(c)];
                    
                    spectra = [spectra ;d];
                    
                end
        else
            if contains(c, "YY")
               c = char("YYYY");
            end
            eval(['ndbc_data.ID_' a '.year_' b '.' c '= d ;']);
        end
       end
       if parameter == "swden"
            si = size(d);
            si2= size(freq);
            spectra = reshape(spectra, [si(1),si2(2)]);
            freq = freq';
            spectra = spectra';
            eval(['ndbc_data.ID_' a '.year_' b '.frequency = freq ;']);
            eval(['ndbc_data.ID_' a '.year_' b '.spectrum = spectra ;']);
        end
        freq = [];
        spectra = [];
        eval(['ndbc_data.ID_' a '.year_' b '.time = datetime(ndbc_data.ID_' a '.year_' b...
            '.YYYY,ndbc_data.ID_' a '.year_' b ...
            '.MM,ndbc_data.ID_' a '.year_' b '.DD, ndbc_data.ID_' a '.year_' b '.hh, 0,0);']);
    end
   end
 
else
    for i=1:length(key)
        sha=cell(datapd{string(py.str(key{i}))}.values.shape);
        x=int64(sha{1,1});
        y=int64(sha{1,2});
        df = reshape(double(py.array.array('d',py.numpy.nditer(datapd{string(py.str(key{i}))}.values,...
    pyargs("flags",{"refs_ok"})))),[x,y]);
        cols = cell(py.list(py.numpy.nditer(datapd{string(py.str(key{i}))}.columns,pyargs("flags",{"refs_ok"}))));
        
        for k = 1:length(cols)
            
            a = char(string(py.str(key{i})));
            b = char(string(py.str(cols{k})));
            c = df(:,k);
            
            if parameter == "swden"
                
                if ~contains(b, ".")
                    if contains(b, "YY")
                        
                        b = char("YYYY");
                        
                    end
                    eval(['ndbc_data.year_' a '.' b '= c ;']);
                    
                else
                    freq = [freq str2double(b)];
                    
                    spectra = [spectra ;c];
                    
                end
            else
               if contains(b, "YY")
                   
                   b = char("YYYY");
               end
               eval(['ndbc_data.year_' a '.' b '= c ;']);
            end
            
        end
        if parameter == "swden"
            si = size(c);
            si2= size(freq);
            spectra = reshape(spectra, [si(1),si2(2)]);
            freq = freq';
            spectra = spectra';
            eval(['ndbc_data.year_' a '.frequency = freq ;']);
            eval(['ndbc_data.year_' a '.spectrum = spectra ;']);
        end
        freq = [];
        spectra = [];
        
        eval(['ndbc_data.year_' a '.time = datetime(ndbc_data.year_' a '.YYYY,ndbc_data.year_' a ...
            '.MM,ndbc_data.year_' a '.DD, ndbc_data.year_' a '.hh,0,0);']);
    end
end


