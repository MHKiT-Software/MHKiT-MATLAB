function solution = minimize_slsqp(inputArg1,inputArg2)
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
% x1: array
%     Component 1 data
% x2: array
%     Component 2 data
% bin_size : int
%     Number of data points in each bin
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


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

