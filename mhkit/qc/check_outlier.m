function results = check_outlier(data,bound,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check or outliers using normalized data within a rolling window
%   Upper and lower bounds in standard deviations
%   Data is normalized using (data-mean)/std
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
%             data.time:   1D array of datetimes or posix times
%
%     bound: cell array of floats
%         [lower bound, upper bound] of standard deviations from mean allowed
%         NaN or py.None can be used for either bound 
%
%     key: string (optional)
%         Data column name or translation dictionary key.  If not specified
%         or set to py.None, all columns are used for test.
%
%     window: int (optional)
%         Size of rolling window (in seconds) used to normalize data
%         default = 3600.  If window is set to py.None, data is normalized 
%         using mean and stddev of entire data set (column by column)
%
%     absolute_value: logical (optional)
%         Use the absolute value of the normalized data, default = 1 (True)
%
%     min_failures: int (optional)
%
%         Minimum number of consecutive failures required for reporting,
%         default = 1
%
%     Must set previous arguments to use later optional arguments
%     (i.e. must set key to use min_failures).
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

 py.importlib.import_module('pecos');

  % check to see if a pandas dataframe or not
  if (isa(data,'py.pandas.core.frame.DataFrame')~=1)
    data=qc_data_to_dataframe(data);
  end
  bound = py.list(bound);
  if nargin == 2
    r = struct(py.pecos.monitoring.check_outlier(data,bound));
  elseif nargin == 3
    r = struct(py.pecos.monitoring.check_outlier(data,bound,...
	      varargin{1}));
  elseif nargin == 4
    r = struct(py.pecos.monitoring.check_outlier(data,bound,...
	      varargin{1},varargin{2}));
  elseif nargin == 5
    r = struct(py.pecos.monitoring.check_outlier(data,bound,...
	      varargin{1},varargin{2},varargin{3}));
  elseif nargin == 6
    r = struct(py.pecos.monitoring.check_outlier(data,bound,...
	      varargin{1},varargin{2},varargin{3},varargin{4}));
  else
    ME = MException('MATLAB:qc_outlier','incorrect number of arguments (2 to 7)');
        throw(ME);
  end

  % Convert to qcdata structure
  results.values=double(r.cleaned_data.values);
  results.mask=int64(r.mask.values);

  % Extract time from index, convert to posix then datetime
  ptime = double(py.array.array('d',py.numpy.nditer(r.cleaned_data.index.values)))/1e9;
  results.time=datetime(ptime,'ConvertFrom','posix');

end 
