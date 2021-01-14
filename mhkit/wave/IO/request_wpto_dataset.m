function datast=request_wpto_dataset(data_type,parameter,lat_lon, years,options)

%%%%%%%%%%%%%%%%%%%%
%     Returns data from the WPTO wave hindcast hosted on AWS at the specified 
%       latitude and longitude points. Visit https://registry.opendata.aws/wpto-pds-us-wave/
%       for more information about the dataset and available 
%       locations and years. 
% 
%         Note: To access the WPTO hindcast data, you will need to configure 
%           h5pyd for data access on HSDS. 
%           See the WPTO_hindcast example notebook for more details.  
%     
% Parameters
% ------------
%     data_type : string
%         data set type of interst
%         Options: 'spatial'
%
%     parameter : string or cell array of strings
%         dataset parameter to be downloaded
%         spatial dataset options: 'directionality_coefficient', 'energy_period',
%                 'maximum_energy_direction,'mean_absolute_period',
%                 'mean_zero-crossing_period', 'omni-directional_wave_power',
%                 'peak_period', 'significant_wave_height',
%                 'spectral_width', 'water_depth'  
%
%     lat_lon: cell array or cell array of cell arrays
%         latitude longitude pairs at which to extract data
%
%     years : cell array 
%          Vector of years to access. The years 1979-2010 available. 
%          Examples: {1996} or {2004,2006,2007}
%
%     tree : str | cKDTree (optional)
%          cKDTree or path to .pkl file containing pre-computed tree
%          of lat, lon coordinates, default = py.None
%
%     unscale : bool (optional)
%          Boolean flag to automatically unscale variables on extraction
%          Default = py.True
%
%     str_decode : bool (optional)
%          Boolean flag to decode the bytestring meta data into normal
%          strings. Setting this to False will speed up the meta data read.
%          Default = py.True
%
%     hsds : bool (optional)
%          Boolean flag to use h5pyd to handle .h5 'files' hosted on AWS
%          behind HSDS. Setting to False will indicate to look for files on 
%          local machine, not AWS. Default = py.True
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

arguments 
    data_type
    parameter
    lat_lon
    years
    options.tree = py.None;
    options.unscale = py.True;
    options.str_decode = py.True;
    options.hsds = py.True;
end

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

datap = py.mhkit.wave.io.hindcast.request_wpto_dataset(data_type,py.list(parameter),...
    lat_lon,py.list(years),pyargs("tree",options.tree,"unscale",options.unscale,...
    "str_decode",options.str_decode,"hsds",options.hsds));

datac=cell(datap);
datapd=datac{1}

datamat=datac{2};
matstr=datamat;
meta_ind = cell(py.list(py.numpy.nditer(matstr.index,pyargs("flags",{"refs_ok"}))));
meta_col = cell(py.list(py.numpy.nditer(matstr.columns,pyargs("flags",{"refs_ok"}))));
data_col = cell(py.list(py.numpy.nditer(datapd.columns,pyargs("flags",{"refs_ok"}))));
meta_val = cell(py.list(py.numpy.nditer(matstr.values,pyargs("flags",{"refs_ok"}))));
disp(data_col{1})
%disp(cell(py.list(py.numpy.nditer(datamat.columns.values,pyargs("flags",{"refs_ok"})))))

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
%datast = vals;
% if ~isempty(fieldnames(matstr))
for k = 1:max(size(data_col))
for j = 1:max(size(meta_ind))
    
    for i=1:si(2)
        dat = char(py.str(data_col{k}))
        disp(class(dat))
        num = dat(end)
        num2 = num+1
        dat = dat(1:end-2)
        
        numst = string(num)
        dat = string(strrep(dat,'-','_')) 
        test='location'+string(numst)
        disp(test)
        disp(dat)
        datast.(test).(dat)=vals(:,i);
        %unit=string(matstr.(test));
        %datast.units.(test)=unit;
end
end
end
%disp(datast)
% else
%     datast.spectrum = vals';
%     for i=1:si(2)
%         temp = [temp, double(py.array.array('d',py.numpy.nditer(vv{i})))];
%     end
%     datast.frequency = temp';
% end
% 
%  for i=1:siti(2)
%     datast.time(i)=datetime(string(py.str(ti{i})),'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS')';
%  end
