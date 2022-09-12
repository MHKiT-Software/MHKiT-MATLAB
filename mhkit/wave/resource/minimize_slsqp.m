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
        n1 = n + 1;
        la = max([1, m]);
        
        r   = zeros([m+n+n+2,1]);
        l   = zeros([(n+1)*(n+2)/2,1]);
        x0_ = zeros([n,1]);
        mu  = zeros([la,1]);
        s   = zeros([n+1,1]);
        u   = zeros([n+1,1]);
        v   = zeros([n+1,1]);
        iw = 1 + la + n1*n/2 + 1 + n + n + n + la + n1 + n1 + n1;
        w_  = w(iw:end);
        goto_code = 0;
        
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
            gocode = 200;
        
        elseif mode == 0
            itermx = iter;
            if acc >= 0.0 
                iexact = 0;
            else
                iexact = 1;
            end 
            acc = abs(acc);
            tol = 10*acc;
            iter = 0;
            ireset = 0;
            n1 = n + 1;
            n2 = n1*n/2;
            n3 = n2 + 1;
            s(1:n) = 0.0;  
            mu(1:m) = 0.0; 
        else
            % call functions at current x
            t = f;
            for j = 1:m
                if j <= meq
                    h1 = c(j);
                else
                    h1 = 0.0;
                end
                t = t + mu(j) * max(-c(j),h1);
            end
            h1 = t - t0;
            if iexact+1 == 1
                if h1 <= h3/10.0 || line > 10
                    goto_code = 500;
                end
                alpha = min(max(h3/(2.0*(h3-h1)),alphamin),alphamax);
                goto_code = 300;
            elseif iexact+1 == 2
                goto_code = 400;
            else
                goto_code = 500;
            end
        end

        % reset bfgs matrix

        if goto_code <= 100            
            if ireset > 5
                mode = check_convergence( ...
                    n,f,f0,x,x0,s,h3,tol,tolf,toldf,toldx,0,8);
                return
            else
                l(1:n2) = 0.0; %<-----|
                %l = dcopy(n2,l(1),0,l,1);
                j = 1;
                for i = 1:n
                    l(j) = 1.0;
                    j = j + n1 - i;
                end
            end
        end

        %  main iteration : search direction, steplength, ldl'-update

        if goto_code <= 200
            mode = 9;
            if iter > itermx
                return
            end
            % search direction as solution of qp - subproblem
            u(1:n) = xl(1:n);
            v(1:n) = xu(1:n);
            %u = dcopy(n,xl,1,u,1);
            %v = dcopy(n,xu,1,v,1);
            u = daxpy(n,-1.0,x,1,u,1);
            v = daxpy(n,-1.0,x,1,v,1);
            h4 = 1.0;

            lsq(m,meq,n,n3,la,l,g,a(:,1:n),c,u,v,s,r,w_,mode);

        end

        %



    end

    function lsq(m,meq,n,nl,la,l,g,a,b,xl,xu,x,y,w,mode)
        %   MINIMIZE with respect to X
        %
        %             ||E*X - F||
        %                             
        %   WITH UPPER TRIANGULAR MATRIX E = +D^.5   *L^T 
        %                     AND VECTOR F = -D^-1/2 *L^-1  *G
        %
        %  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE
        %  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS
        % 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L
        %
        %   SUBJECT TO
        %
        %             A(J)*X - B(J) = 0 ,         J=1,...,MEQ,
        %             A(J)*X - B(J) >=0,          J=MEQ+1,...,M,
        %             XL(I) <= X(I) <= XU(I),     I=1,...,N,
        %     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, 
        %     XU. WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), 
        %     XL(N), XU(N). THE WORKING ARRAY W MUST HAVE AT LEAST THE 
        %     FOLLOWING DIMENSION:
        %     DIM(W) =        (3*N+M)*(N+1)                        for LSQ
        %                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI
        %                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI
        %                      with MINEQ = M - MEQ + 2*N
        %     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE.
        %     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR
        %     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION
        %           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS)
        %     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
        %          MODE=1: SUCCESSFUL COMPUTATION
        %               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
        %               3: ITERATION COUNT EXCEEDED BY NNLS
        %               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
        %               5: MATRIX E IS NOT OF FULL RANK
        %               6: MATRIX C IS NOT OF FULL RANK
        %               7: RANK DEFECT IN HFTI
        n1 = n + 1;
        mineq = m - meq;
        m1 = mineq + n + n;

        %  determine whether to solve problem
        %  with inconsistent linerarization (n2=1) or not (n2=0)
        n2 = n1*n/2 + 1;
        if n2 == nl
            n2 = 0;
        else
            n2 = 1;
        end
        n3 = n - n2;

        % recover matrix e and vector f from l and g
        i2 = 1;
        i3 = 1;
        i4 = 1;
        ie = 1;
        i_f = n*n + 1;
        for i = 1:n3
            i1 = n1 - i;
            diagnal = sqrt(l(i2));
            w(i3) = 0.0;
            w(i3:i1) = w(i3);
            w(i3:end) = dcopy(i1-n2,l(i2:end),1,w(i3:end),n);
            w(i3:end) = dscal(i1-n2,diagnal,w(i3:end),n);
            w(i3) = diagnal;
            w(i_f-1+i) = (g(i)-...
                fort_dot(i-1,w(i4:end),1,w(i_f:end),1))/diagnal;
            i2 = i2 + i1 - n2;
            i3 = i3 + n1;
            i4 = i4 + n;
        end
        if n2 == 1
            w(i3) = l(nl);
            w(i4) = 0.0;
            w(i4:end) = dcopy(n3,w(i4:end),0,w(i4:end),1);
            w(i_f-1+n) = 0.0;
        end
        w(i_f:end) = dscal(n, -1.0, w(i_f:end), 1);

        ic = i_f + n;
        id = ic + meq*n;

        if meq > 0
            % recover matrix c from upper part of a
            for i = 1:meq
                w(ic-1+i:end) = dcopy(n,a(i,:),la,w(ic-1+i:end),meq);
            end

            %  recover vector d from upper part of b

            w(id:end) = dcopy(meq,b,1,w(id:end),1);
            w(id:end) = dscal(meq,-1.0,w(id:end),1);
        end

        ig = id + meq;

        if mineq > 0
            % recover matrix g from lower part of a
            dummy_a = a(:);
            for i = 1:mineq
                w(ig-1+i:end) = ...
                    dcopy(n,dummy_a(meq+i:end),la,w(ig-1+i:end),m1);
            end
        end

        % augment matrix g by +i and -i
        ip = ig + mineq;
        for i = 1:n
            w(ip-1+i) = 0.0;
            w(ip-1+i:end) = dcopy(n,w(ip-1+i:end),0,w(ip-1+i:end),m1);
        end
        w(ip) = 1.0;
        w(ip:end) = dcopy(n,w(ip),0,w(ip:end),m1+1);

        im = ip + n;
        for i = 1:n
           w(im-1+i) = 0.0;
           w(im-1+i:end) = dcopy(n,w(im-1+i:end),0,w(im-1+i:end),m1);
        end 
        w(im) = -1;
        w(im:end) = dcopy(n,w(im:end),0,w(im:end),m1+1);
  
        ih = ig + m1*n;
  
        if  mineq > 0 
           % recover h from lower part of b
           w(ih:end) = dcopy(mineq,b(meq+1:end),1,w(ih:end),1);
           w(ih:end) = dscal(mineq,-1,w(ih:end),1);
        end
  
        % augment vector h by xl and xu
  
        il = ih + mineq;
        w(il:end) = dcopy(n,xl,1,w(il:end),1);
        iu = il + n;
        w(iu:end) = dcopy(n,xu,1,w(iu:end),1);
        w(iu:end) = dscal(n,-1,w(iu:end),1);
  
        iw = iu + n;

        [w(ic:ic+(meq*n)-1), w(id:id+meq-1), w(ie:ie+(n*n)-1),...
         w(i_f:i_f+n-1), w(ig:ig+(m1*n)-1), w(ih:ih+m1-1), w(iw:end),...
         x, xnorm, mode] = ...
            lsei(w(ic:ic+(meq*n)-1),... % c
                 w(id:id+meq-1), ...    % d
                 w(ie:ie+(n*n)-1),...   % e
                 w(i_f:i_f+n-1),...     % f
                 w(ig:ig+(m1*n)-1),...  % g
                 w(ih:ih+m1-1),...      % h
                 max(1,meq),meq,n, n, m1,m1,n,w(iw:end),x);
        %           lc      mc, le,me,lg,mg,n, w,       x
  
    end
  
    function [c,d,e,f,g,h,w,x,xnrm,mode] = ...
            lsei(c,d,e,f,g,h,lc,mc,le,me,lg,mg,n,w,x)
        % dim(c) :   formal (lc,n),    actual (mc,n)
        % dim(d) :   formal (lc  ),    actual (mc  )
        % dim(e) :   formal (le,n),    actual (me,n)
        % dim(f) :   formal (le  ),    actual (me  )
        % dim(g) :   formal (lg,n),    actual (mg,n)
        % dim(h) :   formal (lg  ),    actual (mg  )
        % dim(x) :   formal (n   ),    actual (n   )
        % dim(w) :   2*mc+me+(me+mg)*(n-mc)  for lsei
        %          +(n-mc+1)*(mg+2)+2*mg     for lsi
        % dim(jw):   max(mg,l)

        %c,d,e,f,g,h        
%         c = reshape(w(ic:ic+(mc*n)-1),[mc,n]);
%         d = w(id:id+mc-1);
%         e = reshape(w(ie:ie+(me*n)-1),[me,n]);
%         f = w(i_f:i_f+me-1);
%         g = reshape(w(ig:ig+(mg*n)-1),[mg,n]);
%         h = w(ih:ih+mg-1);

%         iw = ih + mineq + 2*n;
%         w_ = w(iw:end);

        mode = 2;

        if mc <= n
            l = n - mc;
            mc1 = mc + 1;
            iw = (l+1)*(mg+2) + 2*mg + mc;
            ie = iw + mc + 1;
            i_f = ie + me*l;
            ig = i_f + me;
    
            % triangularize c and apply factors to e and g

            for i = 1:mc
                j = min(i+1,lc);
                [c(i,1), w(iw+i), c(j,1)] = ...
                h12(1,    i,  i+1,n,c(i,1),lc,w(iw+i),c(j,1),lc, 1, mc-i);
                [c(i,1), w(iw+i), e] = ...
                h12(2,i,i+1,n,c(i,1),lc,w(iw+i),e,le,1,me);
                [c(i,1), w(iw+i), g] = ...
                h12(2,i,i+1,n,c(i,1),lc,w(iw+i),g,lg,1,mg);
            end

            % solve c*x=d and modify f
            mode = 6;
            for i = 1:mc
                if ( abs(c(i,i))<epmach ) 
                    return
                end
                x(i) = (d(i)-fort_dot(i-1,c(i,1),lc,x,1))/c(i,i);
            end

            mode = 1;
            w(mc1) = 0.0;
            w(mc1:end) = dcopy(mg,w(mc1),0,w(mc1:end),1);

            if mc ~= n
                for i = 1:me
                    w(i_f-1+i) = f(i) - fort_dot(mc,e(i:end),le,x,1);
                end

                % store transformed e & g
                for i = 1:me
                    w(ie-1+i:end) = ...
                        dcopy(l,e(i:end),le,w(ie-1+i:end),me);
                end
                for i = 1:mg
                    w(ig-1+i:end) = ...
                        dcopy(l,g(i:end),lg,w(ig-1+i:end),mg);
                end

                if mg > 0
                    % modify h and solve inequality constrained ls problem
                    for i = 1:mg
                        h(i) = h(i) - fort_dot(mc,g(i:end),lg,x,1);
                    end
                    [w(ie:ie+(me*n)-1),w(i_f:i_f+me-1),w(ig:ig+(mg*n)-1),...
                        h, w(mc1:end), x, xnrm, mode] = ...
                        lsi(w(ie:ie+(me*n)-1),...   e
                            w(i_f:i_f+me-1),...     f
                            w(ig:ig+(mg*n)-1),...   g
                            h,me,me,mg,mg,l,w(mc1:end));
                          % h,le,me,lg,mg,n,w
                    if mc == 0
                        return
                    end
                    t = dnrm2(mc,x,1);
                    xnrm = sqrt(xnrm*xnrm+t*t);
                    if mode ~= 1
                        return
                    end
                else
                    % solve ls without inequality constraints
                    mode = 7;
                    k = max(le,n);
                    t = sqrt(epmach);
                    %call hfti(w(ie),me,me,l,w(if),k,1,t,krank,dum,w,w(l+1))
                    xnrm = dum(1);
                    x(mc1:end) = dcopy(l,w(i_f:end),1,x(mc1:end),1);
                    if krank ~= l
                        return
                    end
                    mode = 1;
                end


            end
        end
    end

    function [e,f,g,h,w,x,xnorm,mode] = lsi(e,f,g,h,le,me,lg,mg,n,w)
        % qr-factors of e and application to f
        e = reshape(e,[le,n]);
        g = reshape(g,[lg,n]);
        for i = 1:n
            j = min(i+1,n);
            [e(i:end),t,~] = h12(1,i,i+1,me,e(i:end),1,t,e(j:end),1,le,n-i);
            [e(i:end),t,f] = h12(2,i,i+1,me,e(i:end),1,t,f,1,1,1);
        end

        % transform g and h to get least distance problem
        mode = 5;
        for i = 1:mg
            for j = 1:n
                if abs(e(j,j)) < 2.22e-16 
                    return
                end
                g(i,j) = (g(i,j)-fort_dot(j-1,g(i:end),lg,e(j:end),1))/e(j,j);
            end
            h(i) = h(i) - fort_dot(n, g(i:end),lg,f,1);
        end

        % solve least distance problem
        [x,xnorm,w,mode] = ldp(g,lg,mg,n,h,w);

    end

    function [x,xnorm,w,mode] = ldp(g,mg,m,n,h,w)
        % Least distance programming routine.
        if n <= 0
            mode = 2;
        else
            % state dual problem
            mode = 1;
            x = 0.0;
            xnorm = 0.0;
            if m ~= 0
                iw = 0;
                for j =1:m
                    for i = 1:n
                        iw = iw + 1;
                        w(iw) = g(j,i);
                    end
                    iw = iw + 1;
                    w(iw) = h(j);
                end
                jf = iw + 1;
                for i = 1:n
                    iw = iw+1;
                    w(iw) = 0.0;
                end
                w(iw+1) = 1.0;
                n1 = n + 1;
                iz = iw + 2;
                iy = iz + n1;
                iwdual = iy + m;
                % solve dual problem
                [w(1:n1*m), w(jf:jf+n1-1), w(iy:iy+m-1), ...
                    rnorm, w(iwdual:iwdual+m-1), w(iz:iz+n1-1), mode] = ...
                    nnls(w(1:n1*m),n1,n1,m,...
                    w(jf:jf+n1-1),w(iwdual:iwdual+m-1),w(iz:iz+n1-1));
                if mode == 1
                    mode = 4;
                    if rnorm > 0.0
                        % compute solution of primal problem
                        fac = 1 - fort_dot(m,h,1,w(iy:end),1);
                        if fac > 0
                            mode = 1;
                            fac = 1/fac;
                            for j = 1:n
                                x(j) = fac*fort_dot(...
                                    m,g(j:end),1,w(iy:end),1);
                            end
                            xnorm = dnrm2(n,x,1);
                            % compute lagrange multipliers for primal 
                            % problem
                            w(1) = 0.0;
                            w = dcopy(m,w(1),0,w,1);
                            w = daxpy(m,fac,w(iy),1,w,1);
                        end
                    end
                end
            end
        end
    end

    function [a, b, x, rnorm, w, zz, mode] = nnls(a,mda,m,n,b,w,zz)
        % Nonnegative least squares algorithm.
        %  Given an m by n matrix, {A}, and an m-vector, {b},
        %  compute an n-vector, {x}, that solves the least squares problem:
        %
        %  {A}{x} = {b} subject to {x} >= 0 

        factor = 0.01;
        goto_code = 100;
        mode = 1;
        a = reshape(a,[mda,n]);

        if m <= 0 || n <= 0
            mode = 2;
            return
        end
        iter = 0;
        itmax = 3*n;

        % initialize the arrays index(1:n) and x(1:n)
        x = 0.0;
        index = 1:n;
        iz2 = n;
        iz1 = 1;
        nsetp = 0;
        npp1 = 1;

        while goto_code > 0
            % ******  main loop begins here  ******
            % quit if all coefficients are already in the solution.
            % or if m cols of a have been triangularized.
            if goto_code <= 100
                if iz1 <= iz2 && nsetp < m
                    % compute components of the dual (negative gradient) 
                    % vector w()
                    for iz = iz1:iz2
                        j = index(iz);
                        sm = 0.0;
                        for l = npp1:m
                            sm = sm + a(l,j) * b(l);
                        end
                        w(j) = sm;
                    end
                    goto_code = 150;
                    if goto_code <= 150
                        % find largest positive w(j)
                        wmax = 0.0;
                        for iz = iz1:iz2
                            j = index(iz);
                            if w(j) > wmax
                                wmax = w(j);
                                izmax = iz;
                            end
                        end
                        % if wmax <= 0. go to termination. this indicates 
                        % satisfaction of the kuhn-tucker conditions.
                        if wmax > 0.0
                            iz = izmax;
                            j = index(iz);
    
                            % the sign of w(j) is ok for j to be moved to 
                            % set p begin the transformation and check new 
                            % diagonal element to avoid near linear 
                            % dependence.
    
                            asave = a(npp1,j);
                            [a(j:end),up,~] = ...
                                h12(1,npp1,npp1+1,m,a(j:end),1,0,0,1,1,0);
                            %                       u         up,c  
                            unorm = 0.0;
                            if nsetp ~= 0
                                for l = 1:nsetp
                                    unorm = unorm + a(l,j)^2;
                                end
                            end
                            unorm = sqrt(unorm);
                            if (unorm+abs(a(npp1,j))*factor) - unorm > 0
                                % col j is sufficiently independent.  
                                % copy b into zz, update zz and solve for 
                                % ztest ( = proposed new value for x(j) ).
                                for l = 1:m
                                    zz(l) = b(l);
                                end
                                [a(j:end),up,zz] = h12(2,npp1,npp1+1,m,...
                                    a(j:end),1,up,zz,1,1,1);
                                ztest = zz(npp1)/a(npp1,j);

                                % Check if ztest is positive
                                if ztest > 0
                                    % the index j=index(iz) has been 
                                    % selected to be moved from set z to 
                                    % set p. update b, update indices, 
                                    % apply householder transformations to 
                                    % cols in new set z, zero subdiagonal 
                                    % elts in col j, set w(j)=0.
                                    for l = 1:m
                                        b(l) = zz(l);
                                    end
                                    index(iz) = index(iz1);
                                    index(iz1) = j;
                                    iz1 = iz1 + 1;
                                    nsetp = npp1;
                                    npp1 = npp1 + 1;

                                    if iz1 <= iz2
                                        for jz = iz1:iz2
                                            jj = index(jz);
                                            [a(j:end),up,~] = ...
                                                h12(2,nsetp,npp1,m,...
                                                a(j:end),1,up,a(jj:end),...
                                                1,mda,1);
                                        end
                                    end

                                    if nsetp ~= m
                                        for l = npp1:m
                                            a(l,j) = 0.0;
                                        end
                                    end

                                    w(j) = 0.0;
                                    % solve the triangular system.
                                    % store the solution temporarily in zz
                                    rtnkey = 1;
                                    goto_code = 300;
                                    continue
                                end % if ztest>0
                            end % if x-unorm > 0

                            % reject j as a candidate to be moved from set 
                            % z to set p. restore a(npp1,j), set w(j)=0., 
                            % and loop back to test dual coeffs again.
                            a(npp1,j) = asave;
                            w(j) = 0.0;
                            goto_code = 150;
                            continue
                        else
                            goto_code = 200;
                            continue
                        end % if wmax > 0
                    end % if code <= 150
                end % if iz1<=iz2 and nsetp<m 
            end % if code <= 100

            % ******  end of main loop  ******
            %
            % come to here for termination.
            % compute the norm of the final residual vector.
            if goto_code <= 200
                sm = 0;
                if npp1 <= m
                    for i = npp1:m
                        sm = sm + b(i)^2;
                    end
                else
                    for j = 1:n
                        w(j) = 0.0;
                    end
                end
                rnorm = sqrt(sm);
                return;
            end % code 200

            % the following block of code is used as an internal subroutine
            % to solve the triangular system, putting the solution in zz()
            if goto_code <= 300
                for l = 1:nsetp
                    ip = nsetp + 1 - l;
                    if l ~= 1
                        for ii = 1:ip
                            zz(ii) = zz(ii) - a(ii,jj)*zz(ip+1);
                        end
                    end
                    jj = index(ip);
                    zz(ip) = zz(ip)/a(ip,jj);
                end

                if rtnkey ~= 1 && rtnkey ~= 2
                    return;
                end
                % ******  secondary loop begins here ******
                %
                % iteration counter.
            
                iter = iter + 1;
                if iter > itmax
                    mode = 3;
                    goto_code = 200;
                    continue
                end

                % see if all new constrained coeffs are feasible.
                % if not compute alpha
                alpha = 2.0;
                for ip = 1:nsetp
                    l = index(ip);
                    if zz(ip) <= 0.0 
                        t = -x(l)/(zz(ip)-x(l));
                        if alpha > t
                            alpha = t;
                            jj = ip;
                        end
                    end
                end
                % if all new constrained coeffs are feasible then 
                % alpha will still = 2.    
                % if so exit from secondary loop to main loop
                if abs(alpha-2) <= 0
                    % ******  end of secondary loop  ******
                    for ip = 1:nsetp
                        i = index(ip);
                        x(i) = zz(ip);
                    end
                    % all new coeffs are positive. 
                    % loop back to beginning.
                    goto_code = 100;
                    continue
                else
                    % otherwise use alpha which will be between 0 and 1
                    % to interpolate between the old x and the new zz.
                    for ip = 1:nsetp
                        l = index(ip);
                        x(l) = x(l) + alpha*(zz(ip)-x(l));
                    end

                    % modify a and b and the index arrays to move 
                    % coefficient i from set p to set z.
                    i = index(jj);
                    while any(x <= 0)
                        x(i) = 0.0;
                        if jj ~= nsetp
                            jj = jj + 1;
                            for j = jj:nsetp
                                ii = index(j);
                                index(j-1) = ii;
                                [cc,ss,a(j-1,ii)] = g1(a(j-1,ii),a(j,ii));
                                a(j,ii) = 0.0;
                                for l = 1:n
                                    if l ~= ii
                                        % apply procedure g2 
                                        temp = a(j-1,l);
                                        a(j-1,l) = cc*temp + ss*a(j,l);
                                        a(j,l) = -ss*temp + cc*a(j,l);
                                    end
                                end
                                % apply procedure g2 (cc,ss,b(j-1),b(j))
                                temp = b(j-1);
                                b(j-1) = cc*temp + ss*b(j);
                                b(j) = -ss*temp + cc*a(j,l);
                            end
                        end

                        npp1 = nsetp;
                        nsetp = nsetp - 1;
                        iz1 = iz1 - 1;
                        index(iz1) = i;

                        % see if the remaining coeffs in set p are feasible
                        % They should be because of the way alpha was 
                        % determined. if any are infeasible it is due to 
                        % round-off error.  any that are nonpositive will 
                        % be set to zero and moved from set p to set z.
                        for jj = 1:nsetp
                            i = index(jj);
                            if x(i) <= 0.0
                                continue
                            end
                        end

                        % copy b( ) into zz( ).  then solve again and 
                        % loop back
                        for i = 1:m
                            zz(i) = b(i);
                        end
                        rtnkey = 2;
                        goto_code = 300;
                        break
                    end % while (350 loop)
                end % abs(alpha-2) <= 0
            end % code 300

        end % outer while
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
        a_eq = zeros([meq,n+1]);
        a_ieq = zeros([mieq,n+1]);
        if ~isempty(cons{1})
            for jj = 1:numel(cons{1})
                a_eq(jj,1:n) = approx_derivative(cons{1}{jj}, x);                
            end
        end
        if ~isempty(cons{2})
            for jj = 1:numel(cons{2})
                a_ieq(jj,1:n) = approx_derivative(cons{2}{jj}, x);                
            end
        end

        % Now combine a_eq and a_ieq into a single a matrix
        if m == 0
            out = zeros([la,n+1]);
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
                    dx_out = dx;
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

    function dy_out = dcopy(n,dx,incx,dy,incy)
        if n <= 0 
            return
        end
        
        if incx == 1 && incy == 1
            m = rem(n, 7);
            if m ~= 0
                for i = 1:m
                    dy(i) = dx(i);
                end
                if n < 7
                    dy_out = dy;
                    return
                end
            end
            mp1 = m + 1;
            for i = mp1:7:n
                dy(i) = dx(i);
                dy(i+1) = dx(i+1);
                dy(i+2) = dx(i+2);
                dy(i+3) = dx(i+3);
                dy(i+4) = dx(i+4);
                dy(i+5) = dx(i+5);
                dy(i+6) = dx(i+6);
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
                dy(iy) = dx(ix);
                ix = ix + incx;
                iy = iy + incy;
            end
        end

        dy_out = dy;
    end

    function mode = check_convergence(n,f,f0,x,x0,s,h3,acc,tolf,toldf,...
        toldx,converged,not_converged)
        
        if h3 < acc
            mode = not_converged;
        else
            % if an are ok then it has converged
            ok = abs(f-f0) < acc;
            if ~ok
                ok = dnrm2(n,s,1) < acc;
            end
            % note that these can be ignored if they are < 0
            if ~ok && tolf >= 0.0
                ok = abs(f) < tolf;
            end
            if ~ok && toldf >= 0.0
                ok = abs(f-f0) < toldf;
            end
            if ~ok && toldx >= 0.0
                xmx0 = x - x0;
                ok = dnrm2(n,xmx0,1) < toldx;
            end

            if ok
                mode = converged;
            else
                mode = not_converged;
            end
        end
    end

    function norm = dnrm2(n,x,incx)
        if n<1 || incx <1
            norm = 0.0;
        elseif n == 1
            norm = abs(x(1));
        else
            scale = 0.0;
            ssq = 1.0;
            % the following loop is equivalent to this call to the lapack
            % auxiliary routine:
            % call dlassq( n, x, incx, scale, ssq )
            for ix = 1:incx:1 + (n-1)*incx
                if abs(x(ix)) > 0.0
                    absxi = abs(x(ix));
                    if scale < absxi
                        ssq = 1.0 + ssq*(scale/absxi)^2;
                        scale = absxi;
                    else
                        ssq = ssq + (absxi/scale)^2;
                    end
                end
            end
            norm = scale*sqrt(ssq);
        end
    end
    
    function [u, up, c] = h12(mode,lpivot,l1,m,u,~,up,c,ice,icv,ncv)        
        if 0 >= lpivot || lpivot >= l1 || l1 > m
            return
        end
        cl = abs(u(lpivot));
        if mode ~= 2
            % construct the transformation.
            for j = l1:m
                cl = max(abs(u(j)),cl);
            end
            if cl <= 0
                return
            end
            clinv = 1.0/cl;
            sm = (u(lpivot)*clinv)^2;
            for j = l1:m
                sm = sm + (u(j)*clinv)^2;
            end
            cl = cl*sqrt(sm);
            if u(lpivot) > 0
                cl = -cl;
            end
            up = u(lpivot) - cl;
            u(lpivot) = cl;
        elseif cl < 0.0
            return
        end

        if ncv > 0
            % apply the transformation i+u*(u**t)/b to c.
            b = up*u(lpivot);
            % b must be nonpositive here
            if b < 0.0
                b = 1/b;
                i2 = 1 - icv + ice*(lpivot-1);
                incr = ice*(l1-lpivot);
                for j = 1:ncv
                    i2 = i2 + icv;
                    i3 = i2 + incr;
                    i4 = i3;
                    sm = c(i2)*up;
                    for i = l1:m
                        sm = sm + c(i3)*u(i);
                        i3 = i3 + ice;
                    end
                    if abs(sm) > 0
                        sm = sm*b;
                        c(i2) =c(i2) + sm*up;
                        for i = l1:m
                            c(i4) = c(i4) + sm*u(i);
                            i4 = i4 + ice;
                        end
                    end
                end
            end
        end
    end

    function [c,s,sig] = g1(a,b)
        % Compute orthogonal rotation matrix.
        if abs(a) > abs(b)
            xr = b/a;
            yr = sqrt(1+xr^2);
            c = sign(a)*(1.0/yr);
            s = c*xr;
            sig = abs(a)*yr;
        else
            if abs(b) > 0
                xr = a/b;
                yr = sqrt(1+xr^2);
                s = sign(b)*(1.0/yr);
                c = s*xr;
                sig = abs(b)*yr;
            else
                sig = 0.0;
                c = 0.0;
                s = 0.0;
            end
        end
    end

%% End Helper Functions

end

