function maep=mean_annual_energy_production_matrix(LM,JM,frequency)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Calculates mean annual energy production (MAEP) from matrix data 
%     along with data frequency in each bin
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
%     frequency: Data frequency for each bin. 
%         Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(Hm0_bins,frequency)
%
%        OR
%
%        structure of form: 
%
%           frequency.values
%
%           frequency.Hm0_bins
%
%           frequency.Te_bins
%         
% Returns
% ---------
%     maep: float
%         Mean annual energy production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit_python_utils');
py.importlib.import_module('mhkit');

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

if (isa(frequency,'py.pandas.core.frame.DataFrame')~=1)
    x=size(frequency.values);
    li=py.list();
    if x(2)>1 
        for i = 1:x(2)
            app=py.list(frequency.values(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
        end
    end
    freqpan=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,py.list(frequency.Hm0_bins),int32(x(2)));
end

maep=double(py.mhkit.wave.performance.mean_annual_energy_production_matrix(LMpan,JMpan,freqpan));