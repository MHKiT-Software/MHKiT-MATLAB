function P=dc_power(voltage,current)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     Calculates the real power from DC voltage and current. 
%     
% Parameters
% ------------
%     voltage: Time series of  measured voltages [V]
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,voltage)
%
%        OR
%
%        structure of form: 
%               voltage.voltage : matrix or vector
%
%               voltage.time : time vector
%
%     current: Time series of current [A]
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,current)
%
%        OR
%
%        structure of form:
%               current.current : matrix or vector
%
%               current.time : time vector
%
% Returns
% ---------
%     P: Structure
%
%
%       P.power [W]
%
%       P.gross: gross power from all lines [W] 
%           
%       P.time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');

py.importlib.import_module('mhkit_python_utils');

if (isa(voltage,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(voltage)==1)
        x=size(voltage.voltage);
        li=py.list();
        if x(2)> 1 
            for i = 1:x(2)
                app=py.list(double(voltage.voltage(:,i)));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
            end
            voltage=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(voltage.time(:,1)),li,int32(x(2)));
        elseif x(2)==1
            voltage=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(voltage.time(:,1)),voltage.voltage,int32(x(2)));
        end
        
    else
        ME = MException('MATLAB:dc_power','voltage needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one or create a structure');
        throw(ME);
    end
end

if (isa(current,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(current)==1)
        x=size(current.current);
        li=py.list();
        if x(2)==3 
            for i = 1:x(2)
                app=py.list(double(current.current(:,i)));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
            end
            current=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(current.time(:,1)),li,int32(x(2)));
        elseif x(2)==1
            current=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(current.time(:,1)),current.current,int32(x(2)));
        end
        
    else
        ME = MException('MATLAB:dc_power','voltage needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one or create a structure');
        throw(ME);
    end
end

p_pd=py.mhkit.wave.performance.dc_power(voltage,current);
vals=double(py.array.array('d',py.numpy.nditer(p_pd.values)));
sha=cell(p_pd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});
vals=reshape(vals,[x,y]);
si=size(vals);
% for i=1:si(2)
%     if i == si(2)
%         P.gross=vals(:,i);
%     else
%         P.power{i}=vals(:,i);
%         %P.power{i}=[P.power{i}];
%     end
%  end
P.power=vals(:,1:end-1);
P.gross=vals(:,end);
P.time=double(py.array.array('d',py.numpy.nditer(p_pd.index))).';

