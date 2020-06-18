function DEL=get_DELs(data,chan_info,options)

%%%%%%%%%%%%%%%%%%%%
%     Calculates the damage equivalent load 
%     
% Parameters
% ------------
%     data: with handles- data.data and data.time  
%         data.data contains a vector or matrix containing time-series 
%         statistics of variables 
%
%     chan_info : vector 
%         channel names to be analyzed and corresponding fatigue slope factor "m"
%         ie. ('TwrBsFxt',4)  
%
%     binNum : int (optional) 
%         number of bins for rainflow counting method (minimum=100)
%         to call: get_DELs(data,chan_info,"binNum",binNum)
%
%     t : double or int (optional)
%         Used to control DEL frequency. Default for 1Hz is 600 seconds for 10min data
%         to call: get_DELs(data,chan_info,"t",t)
%     
% Returns
% ---------
%     DEL: Structure 
%         Damage equivalent load of each specified variable
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data 
    chan_info
    options.binNum {mustBeNumeric} = 100;
    options.t {mustBeNumeric} = 600;
end

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

if (isa(data,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(data)==1)
        dsize=size(data.data);
        li=py.list();
        if dsize(2)>1 
            for i = 1:dsize(2)
                app=py.list(data.data(:,i));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
            end
            datapd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(data.time(:,1),li,int32(dsize(2)));
        elseif x(2)==1
            datapd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(data.time,py.numpy.array(data.data),int32(dsize(2)));
        end
        
    else
        ME = MException('MATLAB:get_DELs','data needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

dfDEL = py.mhkit.loads.get_DELs(datapd,chan_info,pyargs('binNum',options.binNum,'t',options.t);


