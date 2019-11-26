function [lam, x] = power(A, x, cplx, maxit, tol)
% FORMAT [lam, x] = utils.power(A, [x], [maxit], [tol])
% A     - Matrix of handle to a one-argument function
%         This function should return its input size when called without 
%         argument (unless an iinitial value is provided).
% x     - Initial value [rand]
% cplx  - Complex input? [false unless x is complex]
% maxit - Maximum number of iterations [20]
% tol   - Tolerance [1e-4]
% 
% Approximate the leading eigenvalue and eigenvector of a linear operator
% using the Power method.

if isnumeric(A)
    A0  = A;
    dim = [size(A0,2) 1];
    A   = @(x) A0*x;
else
    if nargin < 2 || isempty(x)
        dim = A();
    end
end
if nargin < 5 || isnan(tol)
    tol = 1E-7;
end
if nargin < 4 || isnan(maxit)
    maxit = 100;
end
if nargin < 3 || isnan(cplx)
    cplx = false;
end
if nargin < 2 || isempty(x)
    if cplx
        x = complex(rand(dim),rand(dim));
    else
        x = rand(dim);
    end
end


x0  = x;
x   = A(x);
lam = x0(:)' * x(:);
x = x / norm(x(:));
if lam < 0.0
    x = -x;
end
    
i    = 1;
obj  = eps;
gain = Inf;
while i <= maxit && gain > tol
    
    Ax  = A(x);
    lam = x(:)'*Ax(:);
    nrm = norm(Ax(:));
    x   = Ax / nrm;
    if lam < 0.0
        x = -x;
    end
    
    obj  = [obj nrm+lam];
    gain = abs((obj(end) - obj(end-1))/obj(end-1));
    
    i    = i+1;
    
end

x = Ax / lam;