function bins=bin_statistics(data,bin_against,bin_edges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bins calculated statistics against data signal (or channel) 
%   according to IEC TS 62600-3:2020 ED1.
%    
% Parameters
% ------------
%     data: with handles- data.data and data.time  
%         data.data contains a vector or matrix containing time-series 
%         statistics of variables 
%
%     bin_against: vector
%         Data signal to bin data against (ie. current speed)
%
%     bin_edges: vector
%         Bin edges with consistent step size.
%
%     
% Returns
% ---------
%     bins: structure
%
%
%         bins.averages = means of each bin
%
%         bins.std = standard deviation of each bin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');


fn = fieldnames(data);
si = size(fn);
li=py.list();
li2 = py.list();
 

 for k=1:length(fn)
     
     if ~strcmp(fn{k} , 'time')

            eval(['temp = data.' fn{k} ';' ]); 
            indlen = size(temp);
            app=py.list(temp);
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            li2=py.mhkit_python_utils.pandas_dataframe.lis(li2,fn{k});
        
     end
 end
 ind = 1:1:indlen(1);
datapd=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(ind,li,si(1),pyargs('cols',li2));

stat_py = py.mhkit.loads.bin_statistics(datapd,py.numpy.array(bin_against),py.numpy.array(bin_edges));

averages = double(py.array.array('d',py.numpy.nditer(stat_py{1})));
std = double(py.array.array('d',py.numpy.nditer(stat_py{2})));
sha=cell(stat_py{1}.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
pointer = 1;

for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = averages(pointer:pointer+x-1);
        eval(['bins.averages.' fn{k} '= val ;' ]);
        pointer = pointer +x ;
     end
end
 pointer = 1;
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = std(pointer:pointer+x-1);
        eval(['bins.std.' fn{k} '= val ;' ]);
        pointer = pointer +x; 
     end
 end