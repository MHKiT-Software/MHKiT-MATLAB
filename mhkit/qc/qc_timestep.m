% Wrap QC functions around pecos monitoring check functions
% This all assumes a pandas dataframe with a datetime index
function results = qc_timestep(data, freq)
  results = struct(py.pecos.monitoring.check_timestamp(data,freq))
end 
