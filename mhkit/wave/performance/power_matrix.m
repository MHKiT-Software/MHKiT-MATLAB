function PM=power_matrix(LM,JM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Generates a power matrix from a capture length matrix and wave energy 
%     flux matrix
% 
% Parameters
% ------------
%     LM: Capture Length
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,L)
%
%        OR
%
%        structure of form: 
%
%           LM.values
%
%           LM.stat
%
%           LM.Hm0_bins
%
%           LM.Te_bins
%         
%
%     JM: Wave Energy Flux
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,J)
%
%        OR
%
%        structure of form: 
%
%           JM.values: Wave energy flux matrix
%
%           JM.Hm0_bins
%
%           JM.Te_bins
%         
% Returns
% ---------
%     PM: Structure
%
%
%         PM.values: Power matrix
%
%         PM.stat: statistic of the matrix (i.e. "mean", "max", etc.)
%         (string)
%
%         PM.Hm0_bins
%
%         PM.Te_bins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

% need to add asserts for pandas
if (isa(LM,'py.pandas.core.frame.DataFrame')~=1)
    x=size(LM.values);

    li=py.list();
    if x(2)>1 
        for i = 1:x(2)
            app=py.list(LM.values(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
        end
    end

    LMpan=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,py.list(LM.Hm0_bins),int32(x(2)));
    
end

if (isa(JM,'py.pandas.core.frame.DataFrame')~=1)
    x=size(JM.values);
    li=py.list();
    if x(2)>1 
        for i = 1:x(2)
            app=py.list(JM.values(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
        end
    end
    JMpan=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,py.list(JM.Hm0_bins),int32(x(2)));
end

PMpan=py.mhkit.wave.performance.power_matrix(LMpan,JMpan);

vals=double(py.array.array('d',py.numpy.nditer(PMpan.values)));
sha=cell(PMpan.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[y,x]);
vals=transpose(vals);

PM.values=vals;
PM.stat=LM.stat;
PM.Hm0_bins=double(LM.Hm0_bins);
PM.Te_bins=double(LM.Te_bins);