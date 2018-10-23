function [s,rho,A,ll,llm,llp] = multicoil_infer(varargin)
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

% -------------------------------------------------------------------------
% Helper functions to check input arguments
function ok = isarray(X)
    ok = isnumeric(X) || isa(X, 'file_array');
end
function ok = isboundary(X)
    ok = (isnumeric(X) && isscalar(X) && 0 <= X && X <= 1) || ...
         (ischar(X)    && any(strcmpi(X, {'c','circulant','n','neumann'})));
end
function ok = isrealarray(X)
    function okk = isrealtype(T)
        okk = numel(T) > 7 || strcmpi(T(1:7),'complex');
    end
    if isa(X, 'file_array')
        ok = all(cellfun(@isrealtype, {X.dtype}));
    else
        ok = isreal(X);
    end
end

% -------------------------------------------------------------------------
% Parse input
p = inputParser;
p.FunctionName = 'multicoil_infer';
p.addRequired('CoilImages',                  @isarray);
p.addParameter('SensMaps',      [],          @isarray);
p.addParameter('MeanImage',     [],          @isarray);
p.addParameter('Precision',     1,           @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1,           @isnumeric);
p.addParameter('RegCompFactor', 1,           @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegBoundary',   1,           @isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('CovOptim',      true,        @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Parallel',      0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Tolerance',     1E-3,        @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMax',       100,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('SubIterMax',    10,          @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMin',       1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('SubIterMin',    1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLCond',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrev',        [],          @(X) isnumeric(X));
p.addParameter('Verbose',       0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.parse(varargin{:});
rho           = p.Results.MeanImage;
x             = p.Results.CoilImages;
s             = p.Results.SensMaps;
A             = p.Results.Precision;
reg           = p.Results.RegStructure;
alpha         = p.Results.RegCoilFactor;
gamma         = p.Results.RegCompFactor;
bnd           = p.Results.RegBoundary;
vs            = p.Results.VoxelSize;
optim_cov     = p.Results.CovOptim;
optim         = p.Results.SensOptim;
tol           = p.Results.Tolerance;
itermax       = p.Results.IterMax;
subitermax    = p.Results.SubIterMax;
itermin       = p.Results.IterMin;
subitermin    = p.Results.SubIterMin;
verbose       = p.Results.Verbose;
llm           = p.Results.LLCond;
llp           = p.Results.LLPrior;
ll            = p.Results.LLPrev;
% Nw            = p.Results.Parallel;

% -------------------------------------------------------------------------
% Post-process input
N = size(x,4);

% Precision: default = identity
if numel(A) == 1
    A = A * eye(N);
end
% Reg components: just change reg structure
gamma = padarray(gamma(:)', [0 max(0,2-numel(gamma))], 'replicate', 'post');
% Reg factor: ensure zero sum -> propagate their sum to reg components
alpha = padarray(alpha(:), [max(0,N-numel(alpha)) 0], 'replicate', 'post');
gamma = gamma * sum(alpha);
alpha = alpha/sum(alpha);
% Allocate mean
if isempty(rho)
    rho = zeros(size(x,1),size(x,2),size(x,3));
end
% Allocate sensitivity (TODO: do not allocate phase if not needed?)
if isempty(s)
    s = zeros(size(x,1),size(x,2),size(x,3),N,2);
end

% -------------------------------------------------------------------------
% Hierarchical optimisation

if size(s,5) == 1
    optim_cov   =     [0] * optim_cov;
    optim_mag   =     [1] * optim(1);
    optim_phase =     [0] * optim(2);
    gamma_mag   = 10.^[0] * gamma(1);
    gamma_phase = 10.^[0] * gamma(2);
    tol         = 10.^[0] * tol;
else
%     optim_cov   =     [0 0 1] * optim_cov;
%     optim_mag   =     [0 1 1] * optim(1);
%     optim_phase =     [1 1 1] * optim(2);
%     gamma_mag   = 10.^[0 0 0] * gamma(1);
%     gamma_phase = 10.^[0 0 0] * gamma(2);
%     tol         = 10.^[0 0 0] * tol;
    optim_cov   =     [0] * optim_cov;
    optim_mag   =     [0] * optim(1);
    optim_phase =     [1] * optim(2);
    gamma_mag   = 10.^[0] * gamma(1);
    gamma_phase = 10.^[0] * gamma(2);
    tol         = 10.^[0] * tol;
end

% Remove consecutive repeated combinations
prm = [optim_cov ; optim_mag ; optim_phase ; gamma_mag ; gamma_phase ; tol];
prm(:,all(diff(prm,1,2) == 0)) = [];

optim_cov   = prm(1,:);
optim_mag   = prm(2,:);
optim_phase = prm(3,:);
gamma_mag   = prm(4,:);
gamma_phase = prm(5,:);
tol         = prm(6,:);
stop        = zeros(1,numel(tol));
stop(end)   = 1;

% optim_cov   =     [0 0 0 0 1];
% optim_mag   =     [1 1 1 1 1] * optim(1);
% optim_phase =     [1 1 1 1 1] * optim(2);
% gamma_mag   = 10.^[3 2 1 0 0] * gamma(1);
% gamma_phase = 10.^[3 2 1 0 0] * gamma(2);
% tol         = 10.^[0 0 0 0 0] * tol;
% stop        =     [0 0 0 0 1];

% -------------------------------------------------------------------------
% Time execution
if verbose > 0
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Initial estimates
rho = multicoil_mean_ml(x, s, A, rho, [optim_mag(1) optim_phase(1)]);

if optim(2)
    for i=1:5
        s   = multicoil_init_phase(rho, x, s);
        rho = multicoil_mean_ml(x, s, A, rho, [optim_mag(1) optim_phase(1)]);
    end
end

C   = inv(A);
ldC = spm_matcomp('LogDet', C);
if verbose > 1
    multicoil_plot_mean(rho, C, NaN, vs);
end

if isnan(llp)
    % > Initial log-likelihood (prior term)
    llp = multicoil_ll_prior(s, reg, [gamma_mag(1) gamma_phase(1)], ...
                             alpha, bnd, [optim_mag(1) optim_phase(1)], vs);
end
if isnan(llm)
    % > Initial log-likelihood (cond term)
    llm = multicoil_ll_cond(x,s,rho,A) ...
        - size(x,1)*size(x,2)*size(x,3)*ldC;
end

% -------------------------------------------------------------------------
% Loop
it0    = 0;     % > First iteration used to compute gain denominator
upprm  = true;  % > Did we just update parameters?
for it=1:itermax
    
    % ---------------------------------------------------------------------
    % Update parameters
    if upprm
        upprm = false;
        if verbose > 0
            fprintf('Update parameter values:\n');
            fprintf('- Optim covariance: %d\n', optim_cov(1));
            fprintf('- Optim magnitude:  %d\n', optim_mag(1));
            fprintf('- Optim phase:      %d\n', optim_phase(1));
            fprintf('- Regularisation:   [%g %g %g]\n', reg);
            fprintf('- Reg magnitude:    %g\n', gamma_mag(1));
            fprintf('- Reg phase:        %g\n', gamma_phase(1));
            fprintf('- Tolerance:        %g\n', tol(1));
        end
    end
    if verbose > 0
        fprintf('Iteration %d | Sub-iteration %d\n', it, it-it0);
    end
    
    % ---------------------------------------------------------------------
    % Update mean/covariance (closed-form)
    if optim_cov(1)
        if verbose > 0
            fprintf('> Update Covariance\n');
        end
        rho     = multicoil_mean_ml(x, s, A, rho, [optim_mag(1) optim_phase(1)]);
        [C,A]   = multicoil_cov(rho, x, s);
        ldC = spm_matcomp('LogDet', C);
        if verbose > 1
            multicoil_plot_mean(rho, C, ll, vs);
        end
    end
    
    % ---------------------------------------------------------------------
    % Coil-wise sensitivity update
    for n=1:N % randperm(size(x,4))
        
        % Update mean (ML = closed-form / MAP = Gauss-Newton)
        rho = multicoil_mean_ml(x, s, A, rho, [optim_mag(1) optim_phase(1)]);
     
        if verbose > 1
            multicoil_plot_fit(n, x, s, rho, vs)
        end
        
        % Update sensitivity (Gauss-Newton)
        if verbose > 0
            fprintf('> Update Sensitivity: %2d', n);
        end
        
        [s,llm,llp,ok,ls] = multicoil_sensitivity(...
            rho, x, s, ...
            'Index',         n, ...
            'Precision',     A, ...
            'RegStructure',  reg, ...
            'RegCoilFactor', alpha, ...
            'RegCompFactor', [gamma_mag(1) gamma_phase(1)], ...
            'RegBoundary',   bnd, ...
            'VoxelSize',     vs, ...
            'SensOptim',     [optim_mag(1) optim_phase(1)], ...
            'LLPrior',       llp);
        
        if verbose > 0
            if ok, fprintf(' :D (%d)\n', ls);
            else,  fprintf(' :(\n')
            end
        end
        
        if verbose > 1
            multicoil_plot_fit(n, x, s, rho, vs)
            % , '', sprintf('test/brain/fit_movie_%d.gif', n)
        end
        
    end
    
    % ---------------------------------------------------------------------
    % Center sensitivity fields
    sumsen = sum(bsxfun(@times, numeric(s), reshape(alpha, [1 1 1 N])), 4);
    s(:,:,:,:,:) = bsxfun(@minus, numeric(s), sumsen);
    clear sumsen
    % Should I update llp here?
    
    
    % ---------------------------------------------------------------------
    % Update log-likelihood
    llm = llm - size(x,1)*size(x,2)*size(x,3)*ldC; % Add logDet part
    ll = [ll (llm+llp)];
    if verbose > 1
        multicoil_plot_mean(rho, C, ll, vs);
    end
    
    % ---------------------------------------------------------------------
    % Check gain
    if it > 1
        gain = (ll(end) - ll(end-1))/(max(ll(1:end), [], 'omitnan') - min(ll(1:end), [], 'omitnan'));
        if verbose > 0
            switch sign(gain)
                case  1,   sgn = '+';
                case -1,   sgn = '-';
                case  0,   sgn = '=';
                otherwise, sgn = '';
            end
            fprintf('> Gain: %20.10g (%s)\n', gain, sgn);
        end
        if it >= itermax
            if verbose > 0
                fprintf('Reached maximum number of iterations\n');
            end
            break
        end
        if (it-it0) >= subitermax || ...
           ((it-it0) >= subitermin) && (it >= itermin) && abs(gain) < tol(1)
            if verbose > 0
                fprintf('Converged or reached maximum number of iterations\n');
            end
            if stop(1)
                break
            else
                it0   = it;
                upprm = true;
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % Update parameters
    if upprm
        optim_cov   = optim_cov(2:end);
        optim_mag   = optim_mag(2:end);
        optim_phase = optim_phase(2:end);
        gamma_mag   = gamma_mag(2:end);
        gamma_phase = gamma_phase(2:end);
        tol         = tol(2:end);
        stop        = stop(2:end);
    end
    
end


% -------------------------------------------------------------------------
% Time execution
if verbose > 0
    stop = toc(start);
    fprintf('Processing finished: in %s\n', sec2ydhms(stop));
end

end % < function multicoil_infer