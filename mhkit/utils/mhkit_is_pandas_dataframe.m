function tf = mhkit_is_pandas_dataframe(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if the input is a pandas DataFrame
%
% pandas 3 reports the DataFrame class as pandas.DataFrame, earlier
% versions report pandas.core.frame.DataFrame. Both are accepted.
%
% Parameters
% ------------
%     x : any
%         Input to check
%
% Returns
% ---------
%     tf : logical
%         true if x is a pandas DataFrame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = isa(x, 'py.pandas.DataFrame') || isa(x, 'py.pandas.core.frame.DataFrame');

end
