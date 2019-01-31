function [s,rho,A,ll,llm,llp] = multicoil_infer(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images.
%
% FORMAT [sens,mean,prec,ll] = multicoil_infer(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
%
% KEYWORDS
% --------
% Precision     - Noise precision matrix                    [NaN=estimate]
% RegCoilFactor - Regularisation factor per coil            [1]
% RegCoilComp   - Regularisation factor per Re/Im component [1E5]
% VoxelSize     - Voxel size                                [1 1 1]
% SensOptim     - Optimize real and/or imaginary parts      [true true]
% CovOptim      - Optimize noise covariance                 [false]
% Tolerance     - Convergence threshold                     [1E-3]
% IterMax       - Total maximum number of iterations        [100]
% IterMin       - Total minimum number of iterations        [1]
% Verbose       - (-1=quiet,0=moderate,1=verbose,2=plot)    [0]
%
% ADVANCED
% --------
% SensMaps      - Initial complex sensitivity profiles      (File)Array [Nx Ny Nz Nc (2)]
% MeanImage     - Initial complex mean image                (File)Array [Nx Ny Nz  1 (2)]
% RegStructure  - Regularisation Structure (abs memb bend)  [0 0 1]
% RegBoundary   - Boundary conditions for sensitivities     ['neumann']
% Parallel      - Activate parallelisation                  [false]
% LLCond        - Previous conditional log-likelihood       [NaN=compute]
% LLPrior       - Previous prior log-likelihood             [NaN=compute]
% LLPrev        - Log-likelihood of previous iterations     []
%
% OUTPUT
% ------
% sens - (Log)-Sensitivity maps - (File)Array [Nx Ny Nz Nc (2)]
% mean - Mean image             - (File)Array [Nx Ny Nz Nc (2)]
% prec - Noise precision        -       Array [Nc Nc]
% ll   - Log-likelihood
%
% An output FileArray can be provided by using `MeanImage` as an input. If  
% not provided, the output volume will have the same format as the input  
% coil volume.
%
% Nc = number of coils
% Nx = Phase encode 1 
% Ny = Phase encode 2 /or/ Slice
% Nz = Frequency readout
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% =========================================================================
%
%                       PARSE AND PROCESS ARGUMENTS
%
% =========================================================================

% -------------------------------------------------------------------------
% Helper functions to check input arguments
% -------------------------------------------------------------------------
function ok = isarray(X)
    ok = isnumeric(X) || islogical(X) || isa(X, 'file_array');
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
% -------------------------------------------------------------------------
N = size(varargin{1},4);
p = inputParser;
p.FunctionName = 'multicoil_infer';
p.addRequired('CoilImages',                  @isarray);
p.addParameter('SensMaps',      [],          @isarray);
p.addParameter('MeanImage',     [],          @isarray);
p.addParameter('Precision',     NaN,         @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1/N,         @isnumeric);
p.addParameter('RegCompFactor', 1E5,         @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegBoundary',   1,           @isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('CovOptim',      false,       @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Parallel',      0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('Tolerance',     1E-3,        @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMax',       100,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMin',       1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLCond',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrior',       NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrev',        [],          @(X) isnumeric(X));
p.addParameter('Verbose',       0,           @(X) (isnumeric(X) || islogical(X)) && isscalar(X));
p.addParameter('SamplingMask', [],           @isarray);
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
itermin       = p.Results.IterMin;
verbose       = p.Results.Verbose;
llm           = p.Results.LLCond;
llp           = p.Results.LLPrior;
ll            = p.Results.LLPrev;
Nw            = p.Results.Parallel;
mask          = p.Results.SamplingMask;

% -------------------------------------------------------------------------
% Post-process input
% -------------------------------------------------------------------------
N = size(x,4);

% Precision: default = estimate from magnitude
if isnan(A)
    [~,A] = multicoil_init_cov(x);
end
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
init_rho = false;
if isempty(rho)
    rho = zeros(size(x,1),size(x,2),size(x,3), 'like', x);
    init_rho = true;
end
init_s = false;
if isempty(s)
    s = zeros(size(x,1),size(x,2),size(x,3),N,'like',x);
    init_s = true;
end


% =========================================================================
%
%                           INITIAL ESTIMATES
%
% =========================================================================

% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Compute log-determinant of precision matrix
% -------------------------------------------------------------------------
C   = inv(A);
ldC = spm_matcomp('LogDet', C);

% -------------------------------------------------------------------------
% Initial estimate of the mean
% -------------------------------------------------------------------------
if init_rho
    rho = multicoil_mean_ml(x, s, A, rho, optim);
end

% -------------------------------------------------------------------------
% Initial estimate of the sensitivity phase
% -------------------------------------------------------------------------
% Phase is periodic, and phase wraps should not be penalised. However,
% periodic domains are not handled by spm_field, which deals with
% regularization. To circumvent this issue, I try to never wrap the phase.
% A common problem, though, is that the same target value of pi/-pi might 
% be converged towards from two sides (negative and positive), which
% sometimes leads to bad local minima. The solution I found is to
% initialise the sensitivity phase using the Von Mises mean of pointwise
% differences between the mean and the coil image.
% See: * Bishop's PRML - chapter 2.3.8
%      * https://en.wikipedia.org/wiki/Von_Mises_distribution
% Let's do a few iterations of those to centre the mean image.
if init_s
    for i=1:3
        if verbose > 0
            fprintf('> Init sensitivity\n');
        end
        if optim(1)
            s   = multicoil_init_magnitude(rho, x, s);
        end
        if optim(2)
            s   = multicoil_init_phase(rho, x, s);
        end

        if init_rho
            if verbose > 0
                fprintf('> Update mean\n');
            end
            rho = multicoil_mean_ml(x, s, A, rho, optim);
        end
    end
end

% -------------------------------------------------------------------------
% Log-likelihood
% -------------------------------------------------------------------------

if isnan(llp)
    % > Initial log-likelihood (prior term)
    llp = multicoil_ll_prior(s, reg, gamma, alpha, bnd, optim, vs);
end
if isnan(llm)
    % > Initial log-likelihood (cond term)
    llm = multicoil_ll_cond(x,s,rho,A) ...
        - size(x,1)*size(x,2)*size(x,3)*ldC;
end
% ll = [ll (llm+llp)];

if verbose > 1
    multicoil_plot_mean(rho, C, ll, vs);
end

% -------------------------------------------------------------------------
% Loop
% -------------------------------------------------------------------------
for it=1:itermax
    
    if verbose > -1
        fprintf('Iteration %d\n', it);
    end
    
%     % ---------------------------------------------------------------------
%     % Update mean/covariance (closed-form)
%     % ---------------------------------------------------------------------
%     if optim_cov(1)
%         if verbose > 0
%             fprintf('> Update Covariance\n');
%         end
%         % rho     = multicoil_mean_ml(x, s, A, rho, optim);
%         % rho = multicoil_mean_map(x, s, A, rho, [0 1E-5 0], vs, mask);
%         [C,A]   = multicoil_cov(rho, x, s);
%         ldC = spm_matcomp('LogDet', C);
%         if verbose > 1
%             multicoil_plot_mean(rho, C, ll, vs);
%         end
%     end
    
    % ---------------------------------------------------------------------
    % Coil-wise sensitivity update
    % ---------------------------------------------------------------------
    if isdiag(A)
        % If precision matrix is diagonal, computation of the
        % log-likelihood and derivatives in implemented in a slightly
        % faster way (each sensitivity field only depends on one coil).
        % This means that we need to sum the log-likelihood of each coil to
        % get the complete model log-likelihood.
        llm = 0;
    end
    for n=1:N
        
        % -----------------------------------------------------------------
        % Update mean (Gauss-Newton)
        % -----------------------------------------------------------------
%         if verbose > 0
%             fprintf('> Update mean');
%         end
%         [rho,~,~,ok,ls] = multicoil_mean_map(...
%             x, s, rho, ...
%             'Precision',    A, ...
%             'RegFactor',    0, ...
%             'VoxelSize',    vs, ...
%             'SamplingMask', mask);
%         if verbose > 0
%             if ok, fprintf(' :D (%d)\n', ls);
%             else,  fprintf(' :(\n')
%             end
%         end
%         if verbose > 1
%             multicoil_plot_mean(rho, C, ll, vs);
%         end
        
        if verbose > 1
            % multicoil_plot_mean(rho, C, ll, vs);
            multicoil_plot_fit(n, x, s, rho, mask, vs)
        end
        
        % -----------------------------------------------------------------
        % Update sensitivity (Gauss-Newton)
        % -----------------------------------------------------------------
        if verbose > 0
            fprintf('> Update Sensitivity: %2d', n);
        end
        [s,llm1,llp,ok,ls] = multicoil_sensitivity(...
            rho, x, s, ...
            'Index',         n, ...
            'Precision',     A, ...
            'RegStructure',  reg, ...
            'RegCoilFactor', alpha, ...
            'RegCompFactor', gamma, ...
            'RegBoundary',   bnd, ...
            'VoxelSize',     vs, ...
            'SensOptim',     optim, ...
            'LLPrior',       llp, ...
            'Parallel',      Nw, ...
            'SamplingMask',  mask);
        if isdiag(A)
            llm = llm + llm1;
        else
            llm = llm1;
        end
        if verbose > 0
            if ok, fprintf(' :D (%d)\n', ls);
            else,  fprintf(' :(\n')
            end
        end
        
        if verbose > 1
            multicoil_plot_fit(n, x, s, rho, mask, vs)
            % , '', sprintf('test/brain/fit_movie_%d.gif', n)
        end
        
    end
    
    % Add normalisation term
    llm = llm - sum(mask(:))/numel(mask)*size(x,1)*size(x,2)*size(x,3)*ldC;
    
    % ---------------------------------------------------------------------
    % Center sensitivity fields
    % ---------------------------------------------------------------------
    % We know that, at the optimum, sum{alpha_n*s_n} = 0
    % To converge faster, and to avoid bias towards the initial mean
    % estimate, we enforce this condition at each iteration.
    % Note that the log-likelihood changes slightly and might drop
    % during the first iterations.
    sumsen = zeros(size(s,1),size(s,2),size(s,3),1,size(s,5),'single');
    for n=1:N
        sumsen = sumsen + alpha(n)*single(s(:,:,:,n,:));
    end
    for n=1:N
        s(:,:,:,n,:) = s(:,:,:,n,:) - alpha(n) * sumsen;
    end
    % Should I update the log-likelihood? the prior part is easy to
    % update, but the conditional part also changes, and this is
    % slightly more costly to compute.
    % It does not seem to matter much, so I'm not gonna update it.
    % llpsum = multicoil_ll_prior(sumsen, reg, [gamma_mag(1) gamma_phase(1)], 1, bnd, [optim_mag(1) optim_phase(1)], vs);
    % lp = llp - llpsum;
    clear sumsen
    
    for i=1:5
        if verbose > 0
            fprintf('> Update mean');
        end
        [rho,llm,~,ok,ls] = multicoil_mean_map(...
            x, s, rho, ...
            'Precision',    A, ...
            'RegFactor',    0, ...
            'VoxelSize',    vs, ...
            'SamplingMask', mask);
        if verbose > 0
            if ok, fprintf(' :D (%d)\n', ls);
            else,  fprintf(' :(\n')
            end
        end
        if verbose > 1
            multicoil_plot_mean(rho, C, ll, vs);
        end
   end
   llm = llm - sum(mask(:))/numel(mask)*size(x,1)*size(x,2)*size(x,3)*ldC;
    
    % ---------------------------------------------------------------------
    % Update log-likelihood
    % ---------------------------------------------------------------------
    ll = [ll (llm+llp)];
    if verbose > 1
        multicoil_plot_mean(rho, C, ll, vs);
    end
    
    % ---------------------------------------------------------------------
    % Check gain
    % ---------------------------------------------------------------------
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
            if verbose > -1
                fprintf('Reached maximum number of iterations (%d)\n', it);
            end
            break
        end
        if (it >= itermin) && abs(gain) < tol
            if verbose > -1
                fprintf('Converged (%f < %f)\n', abs(gain), tol);
            end
            break
        end
    end
    
    % TEST: decrease regularisation every 5 iterations
    if mod(it,5) == 0
        gamma = gamma./10;
    end
    
end


% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    stop = toc(start);
    fprintf('Processing finished: in %s\n', sec2ydhms(stop));
end

end % < function multicoil_infer