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
w = zeros([len_w,1]);

% Bounds
xl = NaN([n,1]);
xu = NaN([n,1]);

g = wrapped_grad(x0);
g(end+1) = 0.0;
c = eval_constraints(x0, cons);
a = eval_con_normals(x0, cons, la, n, m, meq, mieq);
fx = obj_func(x0, inp1, inp2);

% Initialize internal SLSQP state variables
iter    = 1;
mode    = 0;
alpha   = 0;
f0      = 0;
gs      = 0;
h1      = 0;
h2      = 0;
h3      = 0;
h4      = 0;
t       = 0;
t0      = 0;
tol     = 0;
iexact  = 0;
incons  = 0;
ireset  = 0;
itermx  = 0;
line    = 0;
n1      = 0;
n2      = 0;
n3      = 0;

while true

    solution = slsqp(m, meq, x0, xl, xu, fx, c, g, a, acc, ...
              iter, mode, len_w, w, alpha, f0, gs, h1, h2, h3, h4, t, t0, ...
              tol, iexact, incons, ireset, itermx, line, n1, n2, n3 );

    if solution.mode == 1
    end

    if solution.mode == -1
    end

    if abs(solution.mode) ~= 1
        break
    end
end

%% SLSQP Routine
    function out = slsqp(m, meq, x, xl, xu, f, c, g, a, acc, ...
              iter, mode, l_w, w, alpha, f0, gs, h1, h2, h3, h4, t, t0, ...
              tol, iexact, incons, ireset, itermx, line, n1, n2, n3 )

        % Parameters
        % ----------
        % m : integer
        %     the total number of constraints
        % meq : integer
        %     the number of equality constraints
        % x : array
        %     stores the current iterate of the n vector x
        % xl : array
        %     an n vector of lower bounds xl to x
        % xu : array
        %     an n vector of upper bounds xu to x
        % f : function handle
        %     the value of the objective function
        % c : array
        %     stores the m vector c of constraints, equality constraints 
        %     (if any) first. Dimension of c must be greater or equal la,
        %     which must be greater or equal max(1,m).
        % g : array
        %     stores the n vector g of partials of the objective function
        % a : array
        %     the m by n matrix a of constraint normals
        % acc : float
        %     controls the final accuracy
        % iter : integer
        %     maximum number of iterations
        % mode : integer
        %     mode controls calculation
        %     **evaluation modes**:
        %    
        %     * ** -1 **: gradient evaluation, (g & a)
        %     * **  0 **: *on entry*: initialization, (f, g, c, a),
        %       *on exit*: required accuracy for solution obtained
        %     * **  1 **: function evaluation, (f & c)
        %    
        %     **failure modes**:
        %    
        %     * ** 2 **: number of equality constraints larger than `n`
        %     * ** 3 **: more than 3*n iterations in [[lsq]] subproblem
        %     * ** 4 **: inequality constraints incompatible
        %     * ** 5 **: singular matrix `e` in [[lsq]] subproblem
        %     * ** 6 **: singular matrix `c` in [[lsq]] subproblem
        %     * ** 7 **: rank-deficient equality constraint subproblem [[hfti]]
        %     * ** 8 **: positive directional derivative for linesearch
        %     * ** 9 **: more than iter iterations in sqp
        %     * ** >=10 **: working space w too small, w should be 
        %                   enlarged to l_w=mode/1000
        % l_w : integer
        %     length of w
        % w : array
        %     w() is a one dimensional working space.
        %     the first m+n+n*n1/2 elements of w must not be
        %     changed between subsequent calls of [[slsqp]].
        %     on return w(1) ... w(m). contain the multipliers
        %     associated with the general constraints, while
        %     w(m+1) ... w(m+n(n+1)/2) store the cholesky factor
        %     l*d*l(t) of the approximate hessian of the
        %     lagrangian columnwise dense as lower triangular
        %     unit matrix l with d in its diagonal and
        %     w(m+n(n+1)/2+n+2 ... w(m+n(n+1)/2+n+2+m+2n)
        %     contain the multipliers associated with all
        %     constraints of the quadratic program finding
        %     the search direction to the solution x*
        % 
        % Returns
        % -------
        % out: structure
        %    Fields:
        %    -----
        %    x      : solution vector
        %    acc    : final accuracy
        %    mode   : error codes
        %    w      : see input description of w

        n = length(x);
        la = max([1, m]);
        
        r   = zeros([m+n+n+2,1]);
        l   = zeros([(n+1)*(n+2)/2,1]);
        x0_ = zeros([n,1]);
        mu  = zeros([la,1]);
        s   = zeros([n+1,1]);
        u   = zeros([n+1,1]);
        v   = zeros([n+1,1]);
        
        if mode < 0
            % call jacobian at current x 
            % update cholesky-factors of hessian matrix by modified bfgs 
            % formula
            for i = 1:n
               u(i,1) = g(i,1) - fort_dot(m,a(:,i),1,r,1) - v(i);
            end
            
            % l'*s
            k = 0;
            for i = 1:n
               h1 = 0;
               k = k + 1;
               for j = i+1:n
                   k = k + 1;
                   h1 = h1 + l(k)*s(j);
               end
               v(i) = s(i) + h1;
            end
            
            % d*l'*s
            k = 1;
            for i = 1:n
                v(i) = l(k)*v(i);
                k = k + n1 - i;
            end 
            
            % l*d*l'*s
            for i = n:-1:1
               h1 = 0;
               k = i;
               for j = 1:i-1
                   h1 = h1 + l(k)*v(j);
                   k = k + n - j;
               end
               v(i) = v(i) + h1;
            end
            
            h1 = fort_dot(n,s,1,u,1);
            h2 = fort_dot(n,s,1,v,1);
            h3 = 0.2*h2;
            if h1 < h3
                h4 = (h2-h3)/(h2-h1);
                h1 = h3;
                u = dscal(n,h4,u,1);
                u = daxpy(n,1-h4,v,1,u,1);
            end

            ldl(n,l,u,+one/h1,v);
            ldl(n,l,v,-one/h2,u);
        
        elseif mode == 0
            something
        else
            % call functions at current x
        end
    end

    function out = ldl(n,a,z,sigma,w)
        ME = MException('MATLAB:wave.resource:minimize_slsqp',...
        "LDL not implemented");
        throwAsCaller(ME);
    end

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
        f0_ = obj_func(x0, inp1, inp2);        
        J_transposed = zeros([length(f0_),length(x0)]);
        h_vecs = diag(h);
        for i = 1:numel(h)
            x = x0 + h_vecs(:,i);
            dx = x(i) - x0(i);
            df = obj_func(x, inp1, inp2) - f0_;
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
        f0_ = fun(x0);        
        J_transposed = zeros([length(f0_),length(x0)]);
        h_vecs = diag(h);
        for i = 1:numel(h)
            x = x0 + h_vecs(:,i);
            dx = x(i) - x0(i);
            df = fun(x) - f0_;
            J_transposed(i) = df / dx;
        end

        out = J_transposed';
    end

    function out = fort_dot(n, dx, incx, dy, incy)        
        dtemp = 0;
        if incx == 1 && incy == 1
            m = rem(n,5);
            if m ~= 0
                for i = 1:m
                    dtempt = dtemp + dx(i)*dy(i);
                end
                if n < 5
                    out = dtemp;
                    return
                end            
            end
            mp1 = m + 1;
            for i = mp1:5:n
                dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + ...
                        dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4);
            end
            out = dtemp;

        else
            ix = 1;
            iy = 1;
            if incx < 0
                ix = (-n+1)*incx + 1;
            end
            if incy < 0
                iy = (-n+1)*incy + 1;
            end
            for i = 1:n
                dtemp = dtemp + dx(ix)*dy(iy);
                ix = ix + incx;
                iy = iy + incy;
            end
            out = dtemp;
        end

    end

    function dx_out = dscal(n, da, dx, incx)        
        if incx == 1
            m = rem(n,5);
            if m ~= 0
                for i = 1:m
                    dx(i) = da*dx(i);
                end
                if n < 5 
                    return; 
                end
            end
            mp1 = m + 1;
            for i = mp1:5:n
                dx(i) = da*dx(i);
                dx(i+1) = da*dx(i+1);
                dx(i+2) = da*dx(i+2);
                dx(i+3) = da*dx(i+3);
                dx(i+4) = da*dx(i+4);
            end

        else
            nincx = n*incx;
            for i = 1:incx:nincx
                dx(i) = da*dx(i);
            end
        end

        dx_out = dx;
    end

    function dy_out = daxpy(n,da,dx,incx,dy,incy)
        if n <= 0 
            return
        end
        if abs(da) <= 0
            return
        end

        if incx == 1 && incy == 1
            m = rem(n, 4);
            if m ~= 0
                for i = 1:m
                    dy(i) = dy(i) + da*dx(i);
                end
                if n < 4
                    dy_out = dy;
                    return
                end
            end
            mp1 = m + 1;
            for i = mp1:4:n
                dy(i) = dy(i) + da*dx(i);
                dy(i+1) = dy(i+1) + da*dx(i+1);
                dy(i+2) = dy(i+2) + da*dx(i+2);
                dy(i+3) = dy(i+3) + da*dx(i+3);
            end
        else
            ix = 1;
            iy = 1;
            if incx < 0
                ix = (-n+1)*incx + 1;
            end
            if incy < 0
                iy = (-n+1)*incy + 1;
            end
            for i = 1:n
                dy(iy) = dy(iy) + da*dx(ix);
                ix = ix + incx;
                iy = iy + incy;
            end
        end

        dy_out = dy;
    end

%% End Helper Functions

end

