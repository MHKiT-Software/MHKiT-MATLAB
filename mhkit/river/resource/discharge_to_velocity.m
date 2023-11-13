function V=discharge_to_velocity(Q,polynomial_coefficients)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates velocity given discharge data and the relationship between
%     discharge and velocity at an individual turbine
%
% Parameters
% ------------
%     Q : Discharge data [m3/s]
%
%         Pandas dataframe indexed by time [datetime or s]:
%
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(timeseries,time,x)
%
%         OR
%
%         structure of form:
%
%            Q.Discharge
%
%            Q.time
%
%     polynomial_coefficients : numpy polynomial
%         List of polynomial coefficients that discribe the relationship between
%         discharge and velocity at an individual turbine
%
% Returns
% ------------
%     V: Structure
%
%
%         V.V: Velocity [m/s]
%
%         V.time: time [datetime or s]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

if (isa(Q,'py.pandas.core.frame.DataFrame')~=1)
    x=size(Q.Discharge);
    li=py.list();
    if x(2)>1
        for i = 1:x(2)
            app=py.list(Q.Discharge(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

        end
    elseif x(2) ==1
        li=Q.Discharge;
    end

    if any(isdatetime(Q.time(1)))
        si=size(Q.time);
        for i=1:si(2)
        Q.time(i)=posixtime(Q.time(i));
        end
    end
    % Q=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,Q.time,int32(x(2)));
    Q = py.mhkit_python_utils.pandas_dataframe.list_to_series(Q.Discharge, Q.time);
end

polynomial_coefficients=py.numpy.poly1d(polynomial_coefficients);

Vdf=py.mhkit.river.resource.discharge_to_velocity(Q,polynomial_coefficients);



xx=cell(Vdf.axes);
v=xx{2};
vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(Vdf.values)));
sha=cell(Vdf.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);

si=size(vals);
 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");

    V.(newname(1))=vals(:,i);

 end
 V.time=double(py.array.array('d',py.numpy.nditer(Vdf.index)));

