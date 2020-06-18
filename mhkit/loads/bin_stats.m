function datast=bin_stats(data,x,bin_edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates bin statistics against "x" according to IEC
%    
% Parameters
% ------------
%     data: with handles- data.data and data.time  
%         data.data contains a vector or matrix containing time-series 
%         statistics of variables 
%
%     x: vector
%         contains array of variable to bin data against (ie. current speed)
%
%     bin_edges: vector
%         contains vector of desired bin edges w/ consistent step size.
%
%     
% Returns
% ---------
%     datast: structure
%
%
%         datast.averages = load means of each bin
%
%         datast.std = additional information related to the binning process 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');


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
        ME = MException('MATLAB:bin_stats','data needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

averages,bstd = py.mhkit.loads.bin_stats(datapd,numpy.array(x),numpy.array(bin_edges));

datast.averages = double(py.array.array('d',py.numpy.nditer(averages)));
datast.std = bstd;