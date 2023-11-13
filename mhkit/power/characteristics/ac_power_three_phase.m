function P=ac_power_three_phase(voltage,current,power_factor,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculates the real power from three phase ac voltage and current.
%
% Parameters
% ------------
%    voltage: Time series of all three measured voltages [V]
%        Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,voltage)
%
%        OR
%
%        structure of form:
%           voltage.voltage : matrix of all three phases
%
%           voltage.time : time vector
%
%    current: Time series of all three measured current [A]
%         Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(time,current)
%
%        OR
%
%        structure of form:
%           current.current : matrix of all three phases
%
%           current.time : time vector
%
%     power_factor : float
%         power factor for the system
%
%     line_to_line: bool (Optional)
%         set true if the given voltage measuremtn is line_to_line
%
% Returns
% ---------
%     P: Structure
%
%
%       P.power [W]
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
        if x(2)==3
            for i = 1:x(2)
                app=py.list(double(voltage.voltage(:,i)));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

            end
            voltage=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(voltage.time(:,1)),li,int32(x(2)));
        elseif x(2)==1
            ME = MException('MATLAB:ac_power_three_phase','Three lines of voltage measurements are required');
            throw(ME);
        end

    else
        ME = MException('MATLAB:ac_power_three_phase','voltage needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one or a structure');
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
            ME = MException('MATLAB:ac_power_three_phase','Three lines of voltage measurements are required');
            throw(ME);
        end

    else
        ME = MException('MATLAB:ac_power_three_phase','voltage needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one or a structure');
        throw(ME);
    end
end

if nargin == 4
    p_pd=py.mhkit.power.characteristics.ac_power_three_phase(voltage,current,power_factor,pyargs('line_to_line',varargin{1}));
elseif nargin ==3
    p_pd=py.mhkit.power.characteristics.ac_power_three_phase(voltage,current,power_factor);
else
    ME = MException('MATLAB:ac_power_three_phase','incorrect number of input arguments');
    throw(ME);
end

P.power=double(py.array.array('d',py.numpy.nditer(p_pd.values)));
P.time=double(py.array.array('d',py.numpy.nditer(p_pd.index)));

