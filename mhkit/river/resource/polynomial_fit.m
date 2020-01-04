function poly=polynomial_fit(x,y,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns a polynomial fit for y given x of order n.
% 
% Parameters
% ----------
%     x : array
%         x data for polynomial fit.
%
%     y : array
%         y data for polynomial fit.
%
%     n : int
%         order of the polynomial fit.
% 
% Returns
% --------
%     poly: structure
%
%
%       poly.coef: polynomial coefficients 
%
%       poly.fit: fit coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
x=py.numpy.array(x);
y=py.numpy.array(y);
n=int32(n);

polyt=py.mhkit.river.resource.polynomial_fit(x,y,n);

polyc=cell(polyt);
coef=polyc{1};
fit=polyc{2};
poly.coef=double(py.array.array('d',py.numpy.nditer(coef.coef)));
poly.fit=fit;


