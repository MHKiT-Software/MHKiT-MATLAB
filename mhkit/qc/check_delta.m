function results = check_delta(data,bound,window,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check for stagnant data and/or abrupt changes in the data using
%   difference between max and min values within a rolling window
%    
% Parameters
% ------------
%
%     data: pandas dataframe or qcdata structure
%          Pandas dataframe indexed by datetime (use 
%          py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(ts,time,x))
%
%          OR
%     
%          qcdata structure of form:
%
%             data.values: 2D array of doubles with arbitrary number of columns
%
%             data.time:   1D array of datetimes or posix times
%
%     bound: cell array of floats
%         [lower bound, upper bound] for min/max delta checking
%         NaN or py.None can be used for either bound 
%
%     window : int or double
%         Size of the rolling window (in seconds) used to compute delta
%
%     key: string (optional)
%         Data column name or translation dictionary key.  If not specified
%         or set to py.None, all columns are used for test.
%         to call: check_delta(data,bound,"key",key)
%
%     direction: string (optional)
%         Options: 'positive', 'negative', or py.None (default)
%           If direction is positive, then only identify positive deltas
%           (the min occurs before the max)
%           If direction is negative, then only identify negative deltas (the max occurs before the min)
%           If direction is py.None, then identify both positive and negative deltas
%         to call: check_delta(data,bound,"direction",direction)
%
%     min_failures: int (optional)
%         Minimum number of consecutive failures required for reporting,
%         default = 1
%         to call: check_delta(data,bound,"min_failures",min_failures)
%
%     
% Returns
% ---------
%     results: qcdata structure of form:
%
%         results.values: array of doubles
%            Same shape as input data.values
%            Elements that failed QC test replaced with NaN 
%
%         results.mask: array of int64
%            Same shape as input data.values
%            Logical mask of QC results (1 = passed, 0 = failed QC test) 
%
%         results.time: array of datetimes
%            Same as input times 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    data
    bound
    window
    options.key = py.None;
    options.direction = py.None;
    options.min_failures = 1;
end

 py.importlib.import_module('pecos');


  % check to see if a pandas dataframe or not
  if (isa(data,'py.pandas.core.frame.DataFrame')~=1)
    data=qc_data_to_dataframe(data);
  end
  bound = py.list(bound);


r = struct(py.pecos.monitoring.check_delta(data,bound,window,...
 	      pyargs("key",options.key,"direction",options.direction,...
          "min_failures",int32(options.min_failures))));

  % Convert to qcdata structure
  results.values=double(r.cleaned_data.values);
  results.mask=int64(r.mask.values);

  % Extract time from index, convert to posix then datetime
  ptime = double(py.array.array('d',py.numpy.nditer(r.cleaned_data.index.values)))/1e9;
  results.time=datetime(ptime,'ConvertFrom','posix');

end 
