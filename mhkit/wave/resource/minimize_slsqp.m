function solution = minimize_slsqp(obj_func, x0, inp1, inp2, constraints)
%MINIMIZE_SLSQP Minimize a scalar function of one or more variables using 
% Sequential Least Squares Programming (SLSQP).
%
%  a nonlinear programming method with quadratic programming subproblems
%  this subroutine solves the general nonlinear programming problem:
%
%  minimize: f(x)
%
%  subject to: 
%      c_j (x) = 0 ,        j = 1,...,meq 
%      c_j (x) >= 0 ,       j = meq+1,...,m 
%      xl_i < x_i < xu_i ,  i = 1,...,n  
%
%  the algorithm implements the method of Han and Powell
%  with BFGS-update of the b-matrix and L1-test function
%  within the steplength algorithm.
%
% Parameters
% ----------
% obj_func : callable
%     The objective function to be minimized.
% 
%         fun(x, *args) -> float
% 
%     where x is an 1-D array with shape (n,) and args
%     are inp1 and inp2, the fixed parameters needed to completely
%     specify the function.
% x0 : ndarray, shape (n,)
%     Initial guess. Array of real elements of size (n,),
%     where n is the number of independent variables.
% inp1 : array
%     argument passed to the objective function and its derivatives
% inp2 : array
%     argument passed to the objective function and its derivatives
% constraints: struct
%     Constraints for SLSQP are defined as a structure.
%     Each field with fields: 
%         type : str
%             Constraint type: 'eq' for equality, 'ineq' for inequality.
%         fun : callable
%             The function defining the constraint.
% 
%     Equality constraint means that the constraint function result is to
%     be zero whereas inequality means that it is to be non-negative.
% 
% Returns
% -------
% solution: structure
%    Fields:
%    -----
%    principal_axes : sign corrected PCA axes
%    shift          : The shift applied to x2
%    x1_fit         : gaussian fit of x1 data
%    mu_param       : fit to _mu_fcn
%    sigma_param    : fit to _sig_fits
%
% Reference
% ----------
%   * Dieter Kraft: "A software package for sequential quadratic 
%     programming", DFVLR-FB 88-28, 1988
%   * https://github.com/jacobwilliams/slsqp/blob/master/src/slsqp_core.f90
%
% License
% ----------
%  Original version copyright 1991: Dieter Kraft, FHM.
%  Released under a BSD license.

maxiter = 100;
acc = 1.0E-6;
epsilon = 1.4901161193847656e-08;
meq = 0;
mieq = 0;
cons_eq = {};
cons_ineq ={};

for field = fieldnames(constraints)'
    if ~isfield(constraints.(field{1}),'type')
        ME = MException('MATLAB:wave.resource:minimize_slsqp',...
        "Constraint %s has no type defined.", field{1});
        throwAsCaller(ME);
    end
    if ~isfield(constraints.(field{1}),'fun')
        ME = MException('MATLAB:wave.resource:minimize_slsqp',...
        "Constraint %s has no function defined.", field{1});
        throwAsCaller(ME);
    end
    if strcmpi(string(constraints.(field{1}).type),"eq")
        meq = meq + 1;
        cons_eq{end+1} = constraints.(field{1}).fun;
    elseif strcmpi(string(constraints.(field{1}).type),"ineq")
        mieq = mieq + 1;
        cons_ineq{end+1} = constraints.(field{1}).fun;
    else
        ME = MException('MATLAB:wave.resource:minimize_slsqp',...
        "Unknown constraint type '%s'.", ...
        string(constraints.(field{1}).type));
        throwAsCaller(ME);
    end
    if ~isa(constraints.(field{1}).fun,'function_handle')
        ME = MException('MATLAB:wave.resource:minimize_slsqp',...
        "The function (fun) for Constraint %s must be a function handle.", ...
        field{1});
        throwAsCaller(ME);
    end
end
cons = {cons_eq, cons_ineq};

% m = The total number of constraints
m = meq + mieq;
% la = The number of constraints, or 1 if there are no constraints
la = max([1, m]);
% n = The number of independent variables
n = length(x0);
% Define the workspaces for SLSQP
n1 = n + 1;
mineq = m - meq + n1 + n1;
len_w = (3*n1+m)*(n1+1)+(n1-meq+1)*(mineq+2) + 2*mineq+(n1+mineq)*(n1-meq) ...
        + 2*meq + n1 + floor(((n+1)*n)/2) + 2*m + 3*n + 3*n1 + 1;
% len_jw = mineq;
% w = zeros(len_w);
% jw = zeros(len_jw);

% Bounds
xl = NaN([n,1]);
xu = NaN([n,1]);

g = wrapped_grad(x0);
g(end+1) = 0.0;
c = eval_constraints(x0, cons);
a = eval_con_normals(x0, cons, la, n, m, meq, mieq);

solution=0;

%% Helper Functions
    function out = eval_constraints(x, cons)
        out = zeros([length(cons{1})+length(cons{2}),1]);
        qq = 1;
        if ~isempty(cons{1})
            for jj = 1:numel(cons{1})
                out(qq) = cons{1}{jj}(x);
                qq = qq + 1;
            end
        end
        if ~isempty(cons{2})
            for jj = 1:numel(cons{2})
                out(qq) = cons{2}{jj}(x);
                qq = qq + 1;
            end
        end
    end 

    function out = eval_con_normals(x, cons, la, n, m, meq, mieq)
        % Compute the normals of the constraints
        a_eq = zeros([meq,n]);
        a_ieq = zeros([mieq,n]);
        if ~isempty(cons{1})
            for jj = 1:numel(cons{1})
                a_eq(jj,:) = approx_derivative(cons{1}{jj}, x);                
            end
        end
        if ~isempty(cons{2})
            for jj = 1:numel(cons{2})
                a_ieq(jj,:) = approx_derivative(cons{2}{jj}, x);                
            end
        end

        % Now combine a_eq and a_ieq into a single a matrix
        if m == 0
            out = zeros([la,n]);
        else
            if meq > 0 && mieq > 0
                out = [a_eq a_ieq];
            elseif meq > 0 && mieq == 0
                out = a_eq;
            else
                out = a_ieq;
            end
        end
    end

    function grad = wrapped_grad(x0)        
        h = ((x0 + 1.4901161193847656e-08) - x0);        
        % Compute finite difference approximation of the derivatives of a
        % vector-valued function.
        % 
        % If a function maps from R^n to R^m, its derivatives form m-by-n 
        % matrix called the Jacobian, where an element (i, j) is a partial 
        % derivative of f[i] with respect to x[j].
        %
        % Only the 2-point method is being implemented 
        f0 = obj_func(x0, inp1, inp2);        
        J_transposed = zeros([length(f0),length(x0)]);
        h_vecs = diag(h);
        for i = 1:numel(h)
            x = x0 + h_vecs(:,i);
            dx = x(i) - x0(i);
            df = obj_func(x, inp1, inp2) - f0;
            J_transposed(i) = df / dx;
        end

        grad = J_transposed';        
    end

    function out = approx_derivative(fun, x0)
        % Compute finite difference approximation of the derivatives of a
        % vector-valued function.
        % 
        % If a function maps from R^n to R^m, its derivatives form m-by-n 
        % matrix called the Jacobian, where an element (i, j) is a partial 
        % derivative of f[i] with respect to x[j].
        %
        % Only the 2-point method is being implemented 
        h = ((x0 + 1.4901161193847656e-08) - x0);
        f0 = fun(x0);        
        J_transposed = zeros([length(f0),length(x0)]);
        h_vecs = diag(h);
        for i = 1:numel(h)
            x = x0 + h_vecs(:,i);
            dx = x(i) - x0(i);
            df = fun(x) - f0;
            J_transposed(i) = df / dx;
        end

        out = J_transposed';
    end
%% End Helper Functions

end

