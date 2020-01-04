function p=velocity_to_power(V,polynomial_coefficients,cut_in,cut_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Calculates power given velocity data and the relationship 
%     between velocity and power from an individual turbine
%     
% Parameters
% ----------
%     V : Velocity [m/s]
%
%          Pandas dataframe indexed by time [datetime or s]  
%
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(timeseries,time,x)
%
%         OR
%
%         structure of form:
%
%           V.V: Velocity [m/s]
%
%           V.time: time [datetime or s]
%
%     polynomial_coefficients : vector
%         vector of polynomial coefficients that discribe the relationship between 
%         velocity and power at an individual turbine
%
%     cut_in: float
%         Velocity values below cut_in are not used to compute P
%
%     cut_out: float
%         Velocity values above cut_out are not used to compute P
%     
% Returns   
% -------
%     p : Structure 
%
%
%        P.P: Power [W] 
%
%        P.time: epoch time [s]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit_python_utils');
py.importlib.import_module('mhkit');

if (isa(V,'py.pandas.core.frame.DataFrame')~=1)
    x=size(V.V);
    li=py.list();
    if x(2)>1 
        for i = 1:x(2)
            app=py.list(V.V(:,i));
            li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
            
        end
    elseif x(2) ==1 
        li=V.V;
    end


    V=py.mhkit_python_utils.pandas_dataframe.timeseries_to_pandas(li,V.time,int32(x(2)));
end

polynomial_coefficients=py.numpy.poly1d(polynomial_coefficients);
cut_in=py.float(cut_in);
cut_out=py.float(cut_out);
Pdf=py.mhkit.river.resource.velocity_to_power(V,polynomial_coefficients,cut_in,cut_out);



xx=cell(Pdf.axes);
v=xx{2};
vv=cell(py.list(py.numpy.nditer(v.values,pyargs("flags",{"refs_ok"}))));

vals=double(py.array.array('d',py.numpy.nditer(Pdf.values)));
sha=cell(Pdf.values.shape);
x=int64(sha{1,1});
y=int64(sha{1,2});

vals=reshape(vals,[x,y]);

si=size(vals);
 for i=1:si(2)
    test=string(py.str(vv{i}));
    newname=split(test,",");
    
    p.(newname(1))=vals(:,i);
    
 end
 p.time=double(py.array.array('d',py.numpy.nditer(Pdf.index)));
