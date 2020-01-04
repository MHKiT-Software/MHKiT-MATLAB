function F=exceedance_probability(Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates the exceedance probability
%    
% Parameters
% ----------
%     Q : Discharge data [m3/s]
%
%         Pandas dataframe indexed by time [datetime or s]  
%
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(timeseries,time,x)
%
%         OR
%
%         structure of form:
%
%           Q.Discharge
%
%           Q.time
%         
% Returns   
% -------
%     F : Structure 
%
%
%         F.F: Exceedance probability [unitless] 
%
%         F.time: time [epoch time (s)]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit_python_utils');
py.importlib.import_module('mhkit');

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
    Q=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,Q.time,int32(x(2)));
end

EPpd=py.mhkit.river.resource.exceedance_probability(Q);

xx=cell(EPpd.axes);
v=xx{2};
vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(EPpd.values)));
sha=cell(EPpd.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);

si=size(vals);
 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");
    
    F.(newname(1))=vals(:,i);
    
 end
 F.time=double(py.array.array('d',py.numpy.nditer(EPpd.index)));




