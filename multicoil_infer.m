function [s,rho,A,ll] = multicoil_infer(x, s, rho, A, prm, vs, itermax, verbose)
% Compute mode estimates (ML, MAP) of the prameters of a probabilistic 
% model of complex multicoil MR images.
%
% FORMAT [s,rho,A,ll] = multicoil_infer(x, s, rho, A, prm, vs, vrb)
%
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% A   -       Array [Nc Nc]           - Noise precision matrix
% prm -       Array [1 3] or [Nc 3]   - Regularisation (/ coil) [a m b]
% vs  -       Array [1 3]             - Voxel size [1 1 1]
% vrb -                               - Verbosity (0=quiet, [1]=print, 2=plot)
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
% An output FileArray can be provided by using `rho` as an input. If not 
% provided, the output volume will have the same format as the input coil 
% volume.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 8
    verbose = 1;
end

% -------------------------------------------------------------------------
% Time execution
if verbose > 0
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Initial mean
rho     = multicoil_mean(x, s, A, rho, [0 1 0], vs);
C       = inv(A);
logdetC = spm_matcomp('LogDet', C);
if verbose > 1
    multicoil_plot_mean(rho, C, NaN, vs);
end

% -------------------------------------------------------------------------
% Loop
do_cov = false; % > Do not try to estimate covariance at first
it0    = 1;     % > First iteration used to compute gain denominator
ll     = NaN;   % > Store log-likelihoods
llp    = 0;     % > Initial log-likelihood (prior term)
llm    = multicoil_ll_cond(x,s,rho,A) - size(x,1)*size(x,2)*size(x,3)*logdetC; % > Initial log-likelihood (cond term)
tol    = 1e-3;  % > Tolerance (if gain lower than this value, stop)
for it=1:itermax
    
    % ---------------------------------------------------------------------
    % Update mean/covariance (closed-form)
    if do_cov
        if verbose > 0
            fprintf('Update Covariance\n');
        end
        rho     = multicoil_mean(x, s, A, rho, [0 1 0], vs);
        [C,A]   = multicoil_cov(rho, x, s);
        logdetC = spm_matcomp('LogDet', C);
        if verbose > 1
            multicoil_plot_mean(rho, C, ll, vs);
        end
    end
    
    % ---------------------------------------------------------------------
    % Coil-wise sensitivity update
    for n=1:size(x,4) % randperm(size(x,4))
        
        % Update mean
        rho = multicoil_mean(x, s, A, rho, [0 1 0], vs);
     
        if verbose > 1
            multicoil_plot_fit(n, x, s, rho, vs)
        end
        
        % Update sensitivity (Gauss-Newton)
        [s,llm,llp] = multicoil_sensitivity(n, rho, x, s, A, prm, vs,llp);
        
        if verbose > 1
            multicoil_plot_fit(n, x, s, rho, vs)
        end
        
    end
    
    % ---------------------------------------------------------------------
    % Update log-likelihood
    llm = llm - size(x,1)*size(x,2)*size(x,3)*logdetC; % Add logDet part
    ll = [ll (llm+llp)];
    if verbose > 1
        multicoil_plot_mean(rho, C, ll, vs);
    end
    
    % ---------------------------------------------------------------------
    % Check gain
    if it > 1
        gain = (ll(end) - ll(end-1))/(max(ll(it0+1:end), [], 'omitnan') - min(ll(it0+1:end), [], 'omitnan'));
        if verbose > 0
            if gain > 0
                sgn = '+';
            elseif gain < 0
                sgn = '-';
            else
                sgn = '=';
            end
            fprintf('Gain: %20.10g (%s)\n', gain, sgn);
        end
        if abs(gain) < tol
            if ~do_cov
                do_cov = true;
                % it0    = it;
                tol    = 0.1 * tol;
            else
                break
            end
        end
    end
end


% -------------------------------------------------------------------------
% Time execution
if verbose > 0
    stop = toc(start);
    fprintf('Processing finished: in %s\n', sec2ydhms(stop));
end