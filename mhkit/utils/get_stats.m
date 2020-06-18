function stats=get_stats(data,freq,varargin)

%%%%%%%%%%%%%%%%%%%%
%     function used to obtain statistics from a dataset 
%     
% Parameters
% ------------
%     data: strucutre  
%         vector of excel datenums to be converted 
%
% Returns
% ---------
%     time: DateTimeIndex 
%         vector of corresponding python datetime values
%
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

fn = fieldnames(data);
si = size(fn);
li=py.list();
li2 = py.list();

 for k=1:length(fn)
     
     if ~strcmp(fn{k} , 'time')

            eval(['temp = data.' fn{k} ]);           
            app=py.list(temp);
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            li2=py.mhkit_python_utils.pandas_dataframe.lis(li2,fn{k});
        %end
     end
 end
 
 datapd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(data.time,li,si(1),pyargs('cols',li2));

stats=1