% Wrap QC functions around pecos monitoring check functions
% This all assumes a pandas dataframe with a datetime index
function results = qc_range(data, range)
  results = struct(py.pecos.monitoring.check_range(data, range))
end 
