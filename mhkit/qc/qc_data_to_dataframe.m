function df = qc_data_to_dataframe(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert qc data structure to pandas dataframe
%    
% Parameters
% ------------
%
%     data: qc data structure
%          structure of form:
%
%             data.values
%             data.time
%     
% Returns
% ---------
%     results: pandas dataframe
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  py.importlib.import_module('pecos');
  py.importlib.import_module('mhkit_python_utils');

  x=size(data.values);
  li=py.list();
  if x(2)>1 
    for i = 1:x(2)
      app=py.list(data.values(:,i));
      li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
    end  
  elseif x(2)==1
    li=data.values;
  end

  % This seems weird, to convert from datetime to posix to python
  % DateTimeIndex, but, whatever, it seems to work, at least for my
  % test case.
  if any(isdatetime(data.time(1)))
    data.time=posixtime(data.time);
  end
  df=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,data.time,int32(x(2)));

  % One final thing to convert index to python DateTimeindex, as required by pecos
  df.index=py.pecos.utils.index_to_datetime(df.index);

end 
