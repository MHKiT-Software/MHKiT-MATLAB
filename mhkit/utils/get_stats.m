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

            eval(['temp = data.' fn{k} ';' ]);           
            app=py.list(temp);
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            li2=py.mhkit_python_utils.pandas_dataframe.lis(li2,fn{k});
        %end
     end
 end
 if any(isdatetime(data.time(1)))
    data.time=posixtime(data.time);
 end
 datapd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(data.time,li,si(1),pyargs('cols',li2));
 datapd.index=py.pecos.utils.index_to_datetime(datapd.index);
 
 if nargin == 3 
     stat_py = py.mhkit.utils.get_stats(datapd,freq,pyargs('period',varargin{1}));
 elseif nargin == 2
      stat_py = py.mhkit.utils.get_stats(datapd,int32(freq)); 
 else
     ME = MException('MATLAB:get_stats','Incorrect number of input arguments');
        throw(ME);
 end
 %names = cell(py.array.tolist(py.numpy.nditer(stat_py{1}.columns.values,pyargs('flags',{'refs_ok'}))))
 mean = double(py.array.array("d",py.numpy.nditer(stat_py{1}.values)));
 max = double(py.array.array("d",py.numpy.nditer(stat_py{2}.values)));
 min = double(py.array.array("d",py.numpy.nditer(stat_py{3}.values)));
 std = double(py.array.array("d",py.numpy.nditer(stat_py{4}.values)));
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = mean(k);
        eval(['stats.mean.' fn{k} '= val ;' ]);
     end
 end
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = max(k);
        eval(['stats.max.' fn{k} '= val ;' ]);
     end
 end
 
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = min(k);
        eval(['stats.min.' fn{k} '= val ;' ]);
     end
 end
 
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = std(k);
        eval(['stats.std.' fn{k} '= val ;' ]);
     end
 end

