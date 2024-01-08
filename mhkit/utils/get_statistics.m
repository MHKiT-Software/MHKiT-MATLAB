function stats=get_statistics(data,freq,options)

%%%%%%%%%%%%%%%%%%%%
%     Calculate mean, max, min and stdev statistics of continuous data for a
%     given statistical window. Default length of statistical window (period) is
%     based on IEC TS 62600-3:2020 ED1. Also allows calculation of statistics for multiple statistical
%     windows of continuous data.
%
% Parameters
% ------------
%     data: strucutre
%         structure of variables to get statistics for with field called
%         time.
%     freq: double or int
%         Sample rate of data [Hz]
%     period: double/int (optional)
%         Statistical window of interest [sec], default = 600
%         to call: get_statistics(data,freq,"period",period)
%     vector_channels: cell array of strings (optional)
%         List of channel names that are to be vector averaged
%         to call: get_statistics(data,freq,"vector_channels",vector_channels)
%
% Returns
% ---------
%     stats: structure
%         Structure with mean, max, min, and stdev of each variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    data
    freq
    options.period = 600;
    options.vector_channels = {};

end

py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

fn = fieldnames(data);
si = size(fn);
li=py.list();
li2 = py.list();

 for k=1:length(fn)

     if ~strcmp(fn{k} , {'time','Timestamp'})

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

stat_py = py.mhkit.utils.get_statistics(datapd,int32(freq),pyargs('period',int32(options.period),'vector_channels',py.list(options.vector_channels)));

mean = double(py.array.array("d",py.numpy.nditer(stat_py{1}.values)));
max = double(py.array.array("d",py.numpy.nditer(stat_py{2}.values)));
min = double(py.array.array("d",py.numpy.nditer(stat_py{3}.values)));
std = double(py.array.array("d",py.numpy.nditer(stat_py{4}.values)));

pointer = 0;

for k=1:length(fn)
     if ~strcmp(fn{k} , {'time','Timestamp'})
        pointer = pointer + 1;
        val1 = mean(pointer);
        val2 = max(pointer);
        val3 = min(pointer);
        val4 = std(pointer);
        eval(['stats.mean.' fn{k} '= val1 ;' ]);
        eval(['stats.max.' fn{k} '= val2 ;' ]);
        eval(['stats.min.' fn{k} '= val3 ;' ]);
        eval(['stats.std.' fn{k} '= val4 ;' ]);
     end
end
