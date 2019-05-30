function [sens,mean,prec,ll,llm,lls] = fit(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images.
%
% FORMAT [sens,mean,prec,ll] = b1m.fit(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
%
% KEYWORDS
% --------
% Precision     - Noise precision matrix                    [NaN=estimate]
% RegCoilFactor - Regularisation factor per coil            [1/Nc]
% RegPartFactor - Regularisation factor per Re/Im part      [1E6 1E6]
% RegDecFactor  - Start with more regularisation            [3] (log10)
% VoxelSize     - Voxel size                                [1]
% SensOptim     - Optimize real and/or imaginary parts      [true true]
% CovOptim      - Optimize noise covariance                 [false]
% Tolerance     - Convergence threshold                     [1E-3]
% IterMax       - Total maximum number of iterations        [15]
% IterMin       - Total minimum number of iterations        [1]
% Verbose       - (-1=quiet,0=moderate,1=verbose,2=plot)    [0]
%
% ADVANCED
% --------
% SamplingMask  - Mask of the k-space sampling scheme       [[]=fully sampled]
% SensMaps      - Initial complex sensitivity profiles      (File)Array [Nx Ny Nz Nc]
% MeanImage     - Initial complex mean image                (File)Array [Nx Ny Nz]
% RegStructure  - Regularisation Structure (abs memb bend)  [0 0 1]
% RegBoundary   - Boundary conditions for sensitivities     ['neumann']
% Parallel      - Number of parallel workers                [0]
% LLCond        - Previous conditional log-likelihood       [NaN=compute]
% LLPrior       - Previous prior log-likelihood             [NaN=compute]
% LLPrev        - Log-likelihood of previous iterations     []
%
% OUTPUT
% ------
% sens - (Log)-Sensitivity maps - (File)Array [Nx Ny Nz Nc]
% mean - Mean image             - (File)Array [Nx Ny Nz Nc]
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
% Parse input
% -------------------------------------------------------------------------
Nc = size(varargin{1},4);
p  = inputParser;
p.FunctionName = 'b1m.fit';
p.addRequired('CoilImages',                  @utils.isarray);
p.addParameter('SensMaps',      [],          @utils.isarray);
p.addParameter('MeanImage',     [],          @utils.isarray);
p.addParameter('Precision',     NaN,         @isnumeric);
p.addParameter('RegStructure',  [0 0 1],     @(X) isnumeric(X) && numel(X) == 3);
p.addParameter('RegCoilFactor', 1/Nc,        @isnumeric);
p.addParameter('RegPartFactor', 1E6,         @(X) isnumeric(X) && numel(X) <= 2);
p.addParameter('RegDecFactor',  3,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('RegBoundary',   1,           @utils.isboundary);
p.addParameter('VoxelSize',     [1 1 1],     @(X) isnumeric(X) && numel(X) <= 3);
p.addParameter('SensOptim',     [true true], @(X) (isnumeric(X) || islogical(X)) && numel(X) == 2);
p.addParameter('CovOptim',      false,       @utils.isboolean);
p.addParameter('Parallel',      0,           @utils.isboolean);
p.addParameter('Tolerance',     1E-3,        @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMax',       15,          @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMin',       1,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLCond',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLSens',        NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrev',        [],          @isnumeric);
p.addParameter('Verbose',       0,           @utils.isboolean);
p.addParameter('SamplingMask',  [],          @utils.isarray);
p.parse(varargin{:});
mean          = p.Results.MeanImage;
coils         = p.Results.CoilImages;
sens          = p.Results.SensMaps;
prec          = p.Results.Precision;
reg           = p.Results.RegStructure;
coilfactor    = p.Results.RegCoilFactor;
partfactor    = p.Results.RegPartFactor;
decfactor     = p.Results.RegDecFactor;
bnd           = p.Results.RegBoundary;
vs            = p.Results.VoxelSize;
optim_cov     = p.Results.CovOptim;
optim         = p.Results.SensOptim;
tol           = p.Results.Tolerance;
itermax       = p.Results.IterMax;
itermin       = p.Results.IterMin;
verbose       = p.Results.Verbose;
llm           = p.Results.LLCond;
lls           = p.Results.LLSens;
ll            = p.Results.LLPrev;
Nw            = p.Results.Parallel;
mask          = p.Results.SamplingMask;

if optim_cov
    warning('Covariance inference is deactivated for now.')
end

% -------------------------------------------------------------------------
% Hard-coded parameters
% -------------------------------------------------------------------------
mean_regfactor = 0;
llr = NaN;

% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Post-process input
% -------------------------------------------------------------------------
% Store a few useful values
Nx   = size(coils,1);
Ny   = size(coils,2);
Nz   = size(coils,3);
Nc   = size(coils,4);
Nvox = Nx*Ny*Nz;
if ~isempty(mask)
    propmask = sum(mask(:))/numel(mask);
else
    propmask = 1;
end
% Precision: default = estimate from magnitude
if numel(prec) == 1
    prec = prec * eye(Nc);
end
% Pad voxel size
vs = padarray(vs(:)', [0 max(0,3-numel(vs))], 'replicate', 'post');
% Reg components: pad reg structure
partfactor = padarray(partfactor(:)', [0 max(0,2-numel(partfactor))], 'replicate', 'post');
% Reg coil factor: ensure zero sum -> propagate their sum to reg components
coilfactor = padarray(coilfactor(:), [max(0,Nc-numel(coilfactor)) 0], 'replicate', 'post');
partfactor = partfactor * sum(coilfactor);
coilfactor = coilfactor/sum(coilfactor);
% Decreasing regularisation
decfactor0 = decfactor;
decfactor  = linspace(decfactor0, 0, itermax)';
decfactor  = 10.^decfactor;
partfactorfinal = partfactor;
% Allocate mean
init_mean = isempty(mean);
if init_mean
    mean = zeros(Nx, Ny, Nz, 'like', coils(1));
end
init_sens = isempty(sens);
if init_sens
    sens = zeros(Nx, Ny, Nz, Nc, 'like', coils(1));
    if isnan(lls)
        lls = 0;
    end
end
% Parallel processing
if isnan(Nw),     Nw = 0; end
if ~isfinite(Nw), Nw = parcluster('local'); Nw = Nw.NumWorkers; end


% =========================================================================
%
%                           INITIAL ESTIMATES
%
% =========================================================================

if verbose > -1
    fprintf('Initialisation\n');
end

% -------------------------------------------------------------------------
% Estimate covariance precision
% -------------------------------------------------------------------------
if isempty(prec) || any(any(isnan(prec)))
    [~,prec] = b1m.init.noise(coils);
    prec = prec / Nvox;
end

% -------------------------------------------------------------------------
% Compute log-determinant of precision matrix
% -------------------------------------------------------------------------
C   = inv(prec);
ldC = utils.matcomp('LogDet', C);

% -------------------------------------------------------------------------
% Initial estimate of the mean
% -------------------------------------------------------------------------
if verbose > 1
    b1m.plot.mean(mean, C, ll, vs);
end
if init_mean
    if verbose > 0, fprintf('> Update mean: '); end
    [mean,llm,llr,ok,ls] = b1m.update.mean(...
        coils, sens, mean,              ...
        'RegFactor',    mean_regfactor, ...
        'llPrior',      llr,            ...
        'Precision',    prec,           ...
        'VoxelSize',    vs,             ...
        'SamplingMask', mask);
    llm = llm - propmask*Nvox*ldC;
    if verbose > 0
        if ok, fprintf(' :D (%d)\n', ls);
        else,  fprintf(' :(\n');
        end
    end
    if verbose > 1
        b1m.plot.mean(mean, C, ll, vs);
    end
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
if verbose > 1
    b1m.plot.mean(mean, C, ll, vs);
end
if init_sens
    lls = 0;
    for i=1:3
        if verbose > 0
            fprintf('> Update sensitivity\n');
        end
        if optim(1)
            sens   = b1m.init.magnitude(mean, coils, sens);
        end
        if optim(2)
            sens   = b1m.init.phase(mean, coils, sens);
        end

        if init_mean
            if verbose > 0, fprintf('> Update mean: '); end
            [mean,llm,llr,ok,ls] = b1m.update.mean(...
                coils, sens, mean,              ...
                'RegFactor',    mean_regfactor, ...
                'llPrior',      llr,            ...
                'Precision',    prec,           ...
                'VoxelSize',    vs,             ...
                'SamplingMask', mask);
            llm = llm - propmask*Nvox*ldC;
            if verbose > 0
                if ok, fprintf(' :D (%d)\n', ls);
                else,  fprintf(' :(\n');
                end
            end
            if verbose > 1
                b1m.plot.mean(mean, C, ll, vs);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Log-likelihood
% -------------------------------------------------------------------------

if isnan(lls)
    % > Initial log-likelihood (prior term)
    lls = b1m.ll.sensitivity(sens, reg, partfactorfinal.*decfactor(1), coilfactor, bnd, optim, vs);
end
if isnan(llm)
    % > Initial log-likelihood (cond term)
    llm = b1m.ll.conditional(coils,sens,mean,prec,mask) - propmask*Nvox*ldC;
end
if isnan(llr)
    llr = 0;
end
ll = [ll (sum(llm)+sum(lls)+llr)];

if verbose > 1
    b1m.plot.mean(mean, C, ll, vs);
end

% -------------------------------------------------------------------------
% Loop
% -------------------------------------------------------------------------
for it=1:itermax
    
    partfactor = partfactorfinal .* decfactor(it);
    if it > 1
        lls = lls * decfactor(it) / decfactor(it-1);
    end
    
    if verbose > -1
        if verbose > 0
            fprintf([repmat('-', [1 78]) '\n']);
        end
        fprintf('Iteration %d\n', it);
        if decfactor0
            fprintf('* reg magnitude = %7.1e\n', partfactor(1));
            fprintf('* reg phase     = %7.1e\n', partfactor(2));
        end
    end
    
    % ---------------------------------------------------------------------
    % Update mean/covariance (closed-form)
    % ---------------------------------------------------------------------
%     if optim_cov
%         if verbose > 0
%             fprintf('> Update Covariance\n');
%         end
%         [C,A] = b1m.update.noise(mean, coils, sens);
%         ldC   = spm_matcomp('LogDet', C);
%         if verbose > 1
%             b1m.plot.mean(rho, C, ll, vs);
%         end
%     end
    
    % ---------------------------------------------------------------------
    % Coil-wise sensitivity update (Gauss-Newton)
    % ---------------------------------------------------------------------
    if verbose > 0, fprintf('> Update Sensitivity:'); end
    if isdiag(prec) && verbose < 3
        % If the precision matrix is diagonal, computation of the
        % log-likelihood and derivatives in implemented in a slightly
        % faster way (each sensitivity field only depends on one coil).
        %
        % * This means that we can distribute the processing of each coil!
        %
        % * This also means that we need to sum the log-likelihood of each 
        %   coil to get the complete model log-likelihood.
        lls0 = lls;
        if numel(lls0) == 1
            lls0 = lls0 * ones(1,Nc);
        end
        lls  = zeros(1,Nc);
        llm  = zeros(1,Nc);
        Ad   = diag(prec);
        parfor(n=1:Nc, Nw)

            [sens(:,:,:,n),llm(n),lls(n),ok,ls] = b1m.update.sensitivity(...
                mean, coils(:,:,:,n), sens(:,:,:,n),...
                'Precision',     Ad(n),             ...
                'RegStructure',  reg,               ...
                'RegCoilFactor', coilfactor(n),     ...
                'RegPartFactor', partfactor,        ...
                'RegBoundary',   bnd,               ...
                'VoxelSize',     vs,                ...
                'SensOptim',     optim,             ...
                'LLPrior',       lls0(n),           ...
                'SamplingMask',  mask);
            if verbose > 0
                if ok, fprintf(' [%2d :D (%d)]', n, ls);
                else,  fprintf(' [%2d :(    ]', n);
                end
            end

        end
        llm = sum(llm);
    else
        if isdiag(prec)
            % We do not parallelise, but still need to sum individual 
            % coil-wise log-likelihoods
            llm = 0;
        end
        for n=1:Nc

            if verbose > 0
                if verbose > 1
                    b1m.plot.fit(n, coils, sens, mean, mask, vs);
                end
                if mod(n,4) == 1, fprintf('\n  | '); end
                fprintf('%2d', n);
            end
            [sens,llm1,lls,ok,ls] = b1m.update.sensitivity(...
                mean, coils, sens,              ...
                'Index',         n,             ...
                'Precision',     prec,          ...
                'RegStructure',  reg,           ...
                'RegCoilFactor', coilfactor,    ...
                'RegPartFactor', partfactor,    ...
                'RegBoundary',   bnd,           ...
                'VoxelSize',     vs,            ...
                'SensOptim',     optim,         ...
                'LLPrior',       sum(lls),      ...
                'SamplingMask',  mask);
            if isdiag(prec)
                llm = llm + llm1;
            else
                llm = llm1;
            end
            if verbose > 0
                if ok, fprintf(' :D (%d) | ', ls);
                else,  fprintf(' :(     | ');
                end
                if verbose > 2
                    b1m.plot.fit(n, coils, sens, mean, mask, vs);
                end
            end

        end
    end
    if verbose > 0, fprintf('\n'); end
    
    % Add normalisation term
    llm = llm - propmask*Nvox*ldC;
    
    % ---------------------------------------------------------------------
    % Center sensitivity fields
    % ---------------------------------------------------------------------
    % We know that, at the optimum, sum{alpha_n*s_n} = 0
    % To converge faster, and to avoid bias towards the initial mean
    % estimate, we enforce this condition at each iteration.
    meansens = zeros(Nx, Ny, Nz, 'like', sens(1));
    for n=1:Nc
        meansens = meansens + coilfactor(n) * single(sens(:,:,:,n));
    end
    for n=1:Nc
        sens(:,:,:,n) = sens(:,:,:,n) - meansens;
    end
    clear meansens
    % Zero-centering the sensitivity changes the conditional and prior
    % terms of the log-likelihood.
    % However, in practice, the impact is small, and updating the 
    % log-likelihood is costly. Therefore, we keep it as is.
    % This means that the log-likelihood changes slightly and might drop
    % during the first iterations.
    
    % ---------------------------------------------------------------------
    % Update mean (Gauss-Newton)
    % ---------------------------------------------------------------------
    if verbose > 0, fprintf('> Update mean:\n '); end
    for i=1:3
        [mean,llm,llr,ok,ls] = b1m.update.mean(...
            coils, sens, mean,              ...
            'RegFactor',    mean_regfactor, ...
            'llPrior',      llr,            ...
            'Precision',    prec,           ...
            'VoxelSize',    vs,             ...
            'SamplingMask', mask);
        if verbose > 0
            if ok, fprintf(' :D (%d) |', ls);
            else,  fprintf(' :(     |');
            end
        end
        if verbose > 1
            b1m.plot.mean(mean, C, ll, vs);
        end
   end
   if verbose > 0, fprintf('\n'); end
   llm = llm - propmask*Nvox*ldC;
    
    % ---------------------------------------------------------------------
    % Update log-likelihood
    % ---------------------------------------------------------------------
    ll = [ll (sum(llm)+sum(lls)+llr)];
    if verbose > 1
        b1m.plot.mean(mean, C, ll, vs);
    end
    
    % ---------------------------------------------------------------------
    % Check gain
    % ---------------------------------------------------------------------
    if it > 1
        llmin = min(ll, [], 'omitnan');
        llmax = max(ll, [], 'omitnan');
        gain  = (ll(end) - ll(end-1))/(llmax-llmin);
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
    
end


% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    stop = toc(start);
    fprintf('Processing finished: in %s\n', utils.sec2ydhms(stop));
end

end % < function multicoil_infer