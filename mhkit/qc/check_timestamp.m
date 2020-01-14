function results = check_timestamp(data, freq, varargin)
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
%         Expected start time.  Default: min(data.time)
%
%     expected_end_time: Timestamp (optional)
%         Expected end time.  Default: max(data.time)
%
%     min_failures: int (optional)
%         Minimum number of consecutive failures required for reporting,
%         default = 1
%
%     exact_times: logical (optional)
%         If 1 (True), times are expected to occur at regular intervals
%         (specified by freq) and data is reindexed to match expected frequency 
%         If 0 (False), times only need to occur once or more within each interval
%         (specified by freq) and data is not reindexed
%
%     Must set previous arguments to use later optional arguments
%     (i.e. must set expected_start_time and expected_end_time to use min_failures).
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

 py.importlib.import_module('pecos');

  % check to see if a pandas dataframe or not
  if (isa(data,'py.pandas.core.frame.DataFrame')~=1)
    data=qc_data_to_dataframe(data);
  end

  if nargin == 2
    r = struct(py.pecos.monitoring.check_timestamp(data,freq));
  elseif nargin == 3
    % Unlike the others, this might need a little work to convert 
    % start and end times into something python can use
    r = struct(py.pecos.monitoring.check_timestamp(data,freq,...
	      varargin{1}));
  elseif nargin == 4
    r = struct(py.pecos.monitoring.check_timestamp(data,freq,...
	      varargin{1},varargin{2}));
  elseif nargin == 5
    r = struct(py.pecos.monitoring.check_timestamp(data,freq,...
	      varargin{1},varargin{2},varargin{3}));
  elseif nargin == 6
    r = struct(py.pecos.monitoring.check_timestamp(data,freq,...
	      varargin{1},varargin{2},varargin{3},varargin{4}));
  else
    ME = MException('MATLAB:qc_timestamp','incorrect number of arguments (2 to 6)');
        throw(ME);
  end

  % Convert to qcdata structure
  results.values=double(r.cleaned_data.values);
  results.mask=int64(r.mask.values);

  % Extract time from index, convert to posix then datetime
  ptime = double(py.array.array('d',py.numpy.nditer(r.cleaned_data.index.values)))/1e9;
  results.time=datetime(ptime,'ConvertFrom','posix');

end 
