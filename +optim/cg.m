function [x,nit,rr] = cg(A,b,x,iM,nit,tol,verbose)
% Conjugate gradient solver.
% CG can be used to solve (large) linear systems: A*x = b
%
% FORMAT [x,nit,nrm,tol] = optim.cg(A,b,[x0],[iM],[nit],[tol],[verbose])
% A       - Function handle for the linear operator A*x (left-hand side)
% b       - Target array (right-hand side)
% x0      - Initial guess [0]
% iM      - Function handle for the (inverse) preconditioning matrix [id]
% nit     - Maximum number of iterations [32]
% tol     - Tolerance for early stopping [1E-3]
% verbose - Verbosity level [0]

% Adapted from Mikael's cg_im_solver.

if nargin < 3 || isempty(x)
    x       = zeros(size(b),'single');
end
if nargin < 4 || isempty(iM)
    iM       = @(y) y;
end
if nargin < 5 || isempty(nit) || isnan(nit)
    nit     = 32;
end
if nargin < 6 || isempty(tol) || isnan(tol)
    tol     = 1e-3;
end
if nargin < 7 || isempty(verbose) || isnan(verbose)
    verbose = true;
end

% Initilisation  
%--------------------------------------------------------------------------
bb = sqrt(sum(b(:).*b(:)));         % Norm of b: sqrt(b'*b)
r  = b - A(x); clear b              % Residual: b - A*x
z  = iM(r);                         % Preconditioned residual

rr = sqrt(sum(r(:).*r(:)));         % Norm of r: sqrt(r'*r)
rz = sum(r(:)'*z(:));               % Inner product of r and z
p     = z;                          % Initial conjugate directions p
beta  = 0;                          % Initial step size

if verbose
    fprintf('%g %g\n', bb, sqrt(rr));
end

% Run algorithm
%--------------------------------------------------------------------------
for j=1:nit
    % Calculate conjugate directions P which defines the direction of descent
    %----------------------------------------------------------------------
    p = z + beta*p;
    clear z

    % Finds the step size of the conj. gradient descent
    %----------------------------------------------------------------------
    Ap    = A(p);
    alpha = rz / sum(p(:).*Ap(:));
    
    % Perform conj. gradient descent, obtaining updated X and R, using the 
    % calculated P and alpha
    %----------------------------------------------------------------------
    x = x + alpha * p; 
    r = r - alpha * Ap;
    clear Ap
    
    % Check convergence
    %---------------------------------------------------------------------- 
    rr  = sqrt(sum(r(:).*r(:)))/bb;    
    if verbose
        fprintf('%g\n', rr);
    end
    if rr < tol
        break;
    end
    
    % Update preconditioned residual
    %----------------------------------------------------------------------   
    z   = iM(r);
    
    % Finds the step size for updating P
    %---------------------------------------------------------------------- 
    rz0  = rz;
    rz   = sum(r(:)'*z(:));
    beta = rz / rz0;
end

nit   = j;