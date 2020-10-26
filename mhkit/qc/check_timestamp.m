function results = check_timestamp(data, freq, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Check time series for missing, non-monotonic, and duplicate timestamps
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
%     freq: int
%         Expected time series frequency, in seconds
%
%     expected_start_time: Timestamp (optional)
%         Expected start time in datetime format.  Default: None
%         to call: check_timestamp(data,freq,"expected_start_time",expected_start_time)
%
%     expected_end_time: Timestamp (optional)
%         Expected end time in datetime format.  Default: None
%         to call: check_timestamp(data,freq,"expected_end_time",expected_end_time)
%
%     min_failures: int (optional)
%         Minimum number of consecutive failures required for reporting,
%         default = 1
%         to call: check_timestamp(data,freq,"min_failures",min_failures)
%
%     exact_times: logical (optional)
%         If py.True, times are expected to occur at regular intervals
%         (specified by freq) and data is reindexed to match expected frequency 
%         If py.False, times only need to occur once or more within each interval
%         (specified by freq) and data is not reindexed
%         to call: check_timestamp(data,freq,"exact_times",exact_times)
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
%            Same as input times (possibly reindexed by exact_times)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments 
    data
    freq
    options.expected_start_time = py.None;
    options.expected_end_time = py.None;
    options.min_failures = 1;
    options.exact_times = py.True;
end
 py.importlib.import_module('pecos');
 py.importlib.import_module('pandas');

  % check to see if a pandas dataframe or not
  if (isa(data,'py.pandas.core.frame.DataFrame')~=1)
    data=qc_data_to_dataframe(data);
  end
  
  if options.expected_start_time ~= py.None
      options.expected_start_time = py.pandas.to_datetime(options.expected_start_time);
  end
  
  if options.expected_end_time ~= py.None
      options.expected_end_time = py.pandas.to_datetime(options.expected_end_time);
  end


r = struct(py.pecos.monitoring.check_timestamp(data,freq,...
 	      pyargs("expected_start_time",options.expected_start_time,...
          "expected_end_time",options.expected_end_time,"min_failures",...
          int32(options.min_failures),"exact_times",options.exact_times)));

  % Convert to qcdata structure
  results.values=double(r.cleaned_data.values);
  results.mask=int64(r.mask.values);

  % Extract time from index, convert to posix then datetime
  ptime = double(py.array.array('d',py.numpy.nditer(r.cleaned_data.index.values)))/1e9;
  results.time=datetime(ptime,'ConvertFrom','posix');

end 
