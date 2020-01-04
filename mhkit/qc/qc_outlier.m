% Wrap QC functions around pecos monitoring check functions
% This all assumes a pandas dataframe with a datetime index
function results = qc_corrupt(data,vals)
  results = struct(py.pecos.monitoring.check_corrupt(data,vals))
end 
