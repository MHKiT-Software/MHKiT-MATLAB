function datast=bin_statistics(data,x,bin_edges)

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

stat_py = py.mhkit.loads.bin_statistics(datapd,py.numpy.array(x),py.numpy.array(bin_edges));

averages = double(py.array.array('d',py.numpy.nditer(stat_py{1})));
std = double(py.array.array('d',py.numpy.nditer(stat_py{2})));
sha=cell(stat_py{1}.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
pointer = 1;

for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = averages(pointer:pointer+x-1);
        eval(['datast.averages.' fn{k} '= val ;' ]);
        pointer = pointer +x ;
     end
end
 pointer = 1;
 for k=1:length(fn)
     if ~strcmp(fn{k} , 'time')
        val = std(pointer:pointer+x-1);
        eval(['datast.std.' fn{k} '= val ;' ]);
        pointer = pointer +x; 
     end
 end