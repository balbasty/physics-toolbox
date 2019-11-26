function [sens,meanim,prec,ll,llm,lls,llr] = fit(varargin)
% Compute mode estimates (ML, MAP) of the parameters of a probabilistic 
% model of complex multicoil MR images.
%
% FORMAT [sens,mean,prec,ll] = b1m.fit(coils, ...)
%
% REQUIRED
% --------
% coils - (File)Array [Nx Ny Nz Nch Nct] - Complex coil images
%
% KEYWORDS
% --------
% VoxelSize         - Voxel size                                [1]
% SamplingMask      - Mask of the k-space sampling scheme       [[]=fully sampled]
% Precision         - Noise precision matrix                    [NaN=estimate from fully-sampled]
% SensRegFactor     - Regularisation factor                     [1E2]
% SensRegDecFactor  - Start with more regularisation            [3] (log10)
% MeanRegFactor     - Regularisation factor                     [0]
% MeanRegDecFactor  - Start with more regularisation            [3] (log10)
% IterMax           - Total maximum number of iterations        [15]
% IterDec           - Number of decreasing reg iterations       [10]
% IterMin           - Total minimum number of iterations        [15]
% Tolerance         - Convergence threshold                     [1E-3]
% Verbose           - -1  = quiet
%                     [0] = moderate
%                      1  = verbose
%                      2  = plot mean image and log-likelihood
%                      3  = plot coil-wise fit (slow + disables parallel processing)
%
% ADVANCED
% --------
% SensMaps          - Initial complex sensitivity profiles      (File)Array [Nx Ny Nz Nch]
% MeanImage         - Initial complex mean image                (File)Array [Nx Ny Nz]
% SensOptim         - Optimize sensitivities                    [true]
% MeanOptim         - Optimize mean image                       [true]
% SensLog           - Encode and penalise log-sensitivities     [false]
% CovOptim          - Optimize noise covariance                 [false]
% MeanIter          - Number of MM-Newton updates               [NaN=guess]
% MeanIterInit      - Number of initial MM-Newton updates       [NaN=guess]
% Parallel          - Number of parallel workers                [Inf=all]
% LLCond            - Previous conditional log-likelihood       [NaN=compute]
% LLSens            - Previous sensitivity prior log-likelihood [NaN=compute]
% LLMean            - Previous mean image prior log-likelihood  [NaN=compute]
% LLPrev            - Log-likelihood of previous iterations     [NaN=compute]
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
% Nch = number of coils
% Nx  = Phase encode 1 
% Ny  = Phase encode 2 /or/ Slice
% Nz  = Frequency readout
% Nct = number of contrasts
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
p.addRequired('CoilImages',                      @utils.isarray);
p.addParameter('SensMaps',          [],          @utils.isarray);
p.addParameter('MeanImage',         [],          @utils.isarray);
p.addParameter('Precision',         NaN,         @isnumeric);
p.addParameter('SensRegFactor',     1E4,         @(X) isnumeric(X) && size(X,1) <= Nc && size(X,2) <= 2);
p.addParameter('SensRegDecFactor',  3,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('MeanRegFactor',     0,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('MeanRegDecFactor',  3,           @(X) isnumeric(X) && isscalar(X));
p.addParameter('VoxelSize',         [1 1 1],     @(X) isnumeric(X) && isrow(X) && numel(X) <= 3);
p.addParameter('SensOptim',         true,        @utils.isboolean);
p.addParameter('MeanOptim',         true,        @utils.isboolean);
p.addParameter('SensLog',           false,       @utils.isboolean);
p.addParameter('precOptim',         false,       @utils.isboolean);
p.addParameter('Parallel',          Inf,         @utils.isboolean);
p.addParameter('Tolerance',         1E-3,        @(X) isnumeric(X) && isscalar(X));
p.addParameter('IterMax',           15,          @utils.isintval);
p.addParameter('IterMin',           15,          @utils.isintval);
p.addParameter('IterDec',           10,          @utils.isintval);
p.addParameter('MeanIter',          NaN,         @(X) isscalar(X) && (utils.isintval(X) || isnan(X)));
p.addParameter('MeanIterInit',      NaN,         @(X) isscalar(X) && (utils.isintval(X) || isnan(X)));
p.addParameter('MaskBackground',    false,       @(X) utils.isarray(X) || utils.isboolean(X));
p.addParameter('LLCond',            NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLSens',            NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLMean',            NaN,         @(X) isnumeric(X) && isscalar(X));
p.addParameter('LLPrev',            [],          @isnumeric);
p.addParameter('Verbose',           0,           @utils.isboolean);
p.addParameter('SamplingMask',      [],          @utils.isarray);
p.parse(varargin{:});
meanim          = p.Results.MeanImage;
coils           = p.Results.CoilImages;
sens            = p.Results.SensMaps;
prec            = p.Results.Precision;
sensfactor      = p.Results.SensRegFactor;
sensdecfactor   = p.Results.SensRegDecFactor;
sensoptim       = p.Results.SensOptim;
senslog         = p.Results.SensLog;
meanfactor      = p.Results.MeanRegFactor;
meandecfactor   = p.Results.MeanRegDecFactor;
meanoptim       = p.Results.MeanOptim;
meaniter        = p.Results.MeanIter;
meaniter0       = p.Results.MeanIterInit;
precoptim       = p.Results.precOptim;
vs              = p.Results.VoxelSize;
tol             = p.Results.Tolerance;
itermax         = p.Results.IterMax;
itermin         = p.Results.IterMin;
iterdec         = p.Results.IterDec;
verbose         = p.Results.Verbose;
llm             = p.Results.LLCond;
lls             = p.Results.LLSens;
llr             = p.Results.LLMean;
ll              = p.Results.LLPrev;
Nw              = p.Results.Parallel;
mask            = p.Results.SamplingMask;
maskbg          = p.Results.MaskBackground;

if precoptim
    warning('variance inference is deactivated for now.')
end

% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    fprintf('Processing started\n');
    start = tic;
end

% -------------------------------------------------------------------------
% Store a few useful values
% -------------------------------------------------------------------------
Nx   = size(coils,1);
Ny   = size(coils,2);
Nz   = size(coils,3);
Nc   = size(coils,4);
Nct  = size(coils,5);
Nvox = Nx*Ny*Nz;

% -------------------------------------------------------------------------
% Acceleration-specific parameters
% -------------------------------------------------------------------------
mask = logical(mask);
if isscalar(mask)
    if mask
        mask = logical([]);
    else
        error('The sampling mask is full of zeros')
    end
end
if ~isempty(mask)
    propmask = sum(mask(:))/numel(mask);
else
    propmask = 1;
end
acceleration = 1/propmask;
if acceleration > 1
    if ~isfinite(meaniter)
        meaniter = ceil(acceleration);
    end
    if ~isfinite(meaniter0)
        meaniter0 = 2*ceil(acceleration);
    end
else
    meaniter  = 1;
    meaniter0 = 1;
end

% -------------------------------------------------------------------------
% Precision: default = estimate from magnitude
% -------------------------------------------------------------------------
if numel(prec) == 1
    prec = prec * eye(Nc);
end

% -------------------------------------------------------------------------
% Pad voxel size
% -------------------------------------------------------------------------
vs = utils.pad(vs(:)', [0 3-numel(vs)], 'replicate', 'post');

voxvol = prod(vs);
sensfactor = sensfactor * voxvol;
meanfactor = meanfactor * voxvol;

% -------------------------------------------------------------------------
% Coil factor
% -------------------------------------------------------------------------
sensfactor = utils.pad(double(sensfactor), [Nc-size(sensfactor,1) 2-size(sensfactor,2)], 'replicate', 'post');

% -------------------------------------------------------------------------
% Correct mean image regularisation
% -------------------------------------------------------------------------
meanval = 0;
nbval   = 0;
for n=1:Nc
    coil1   = abs(coils(:,:,:,n));
    coil1   = coil1(:);
    meanval = meanval + sum(coil1, 'omitnan', 'double');
    nbval   = nbval + sum(~isnan(coil1));
    clear coil1
end
meanval    = meanval / nbval;
meanfactor = meanfactor / meanval^2;

% -------------------------------------------------------------------------
% Decreasing regularisation
% -------------------------------------------------------------------------
sensdecfactor0      = sensdecfactor;
sensdecfactor       = linspace(sensdecfactor0, 0, iterdec-1)';
if isempty(sensdecfactor)
    sensdecfactor = 0;
end
sensdecfactor       = 10.^sensdecfactor;
sensfactorfinal     = sensfactor;
sensfactor          = sensdecfactor(1) * sensfactorfinal; 
meandecfactor0      = meandecfactor;
meandecfactor       = linspace(meandecfactor0, 0, iterdec-1)';
if isempty(meandecfactor)
    meandecfactor = 0;
end
meandecfactor       = 10.^meandecfactor;
meanfactorfinal     = meanfactor;
meanfactor          = meandecfactor(1) * meanfactorfinal; 

% -------------------------------------------------------------------------
% Allocate mean
% -------------------------------------------------------------------------
meaninit = isempty(meanim) || (acceleration > 1);
if isempty(meanim)
    meanim = zeros(Nx, Ny, Nz, 1, Nct, 'like', coils(1));
    if isnan(llr)
        llr = 0;
    end
end
if isempty(sens)
    sensinitmag   = true;
    sensinitphase = true;
    sensinit      = true;
    sens = ones(Nx, Ny, Nz, Nc, 'like', coils(1));
    if isnan(lls)
        lls = 0;
    end
else
    sensinitmag   = false;
    sensinitphase = false;
    sensinit      = false;
    if ~isreal(coils) && isreal(sens)
        sensinitphase = true;
        meaninit      = true;
    end
end

% -------------------------------------------------------------------------
% Parallel processing
% -------------------------------------------------------------------------
if isnan(Nw),     Nw = 0; end
if ~isfinite(Nw), Nw = parcluster('local'); Nw = Nw.NumWorkers; end

% -------------------------------------------------------------------------
% Print parameters
% -------------------------------------------------------------------------
if verbose > -1
    if verbose > 0
        fprintf([repmat('-', [1 78]) '\n']);
    end
    fprintf('Parameters\n');
    fprintf(' * Dimensions     | matrix (k1/k2/rd):     [%d %d %d]\n', Nx, Ny, Nz);
    fprintf(' * Dimensions     | coils:                 %d\n', Nc);
    fprintf(' * Dimensions     | contrasts:             %d\n', Nct);
    fprintf(' * Sensitivities  | log-encoding:          %d\n', senslog);
    fprintf(' * Sensitivities  | optimise:              %d\n', sensoptim);
    if sensoptim
        fprintf(' * Sensitivities  | factor:                %g\n', mean(sensfactorfinal(:)));
        fprintf(' * Sensitivities  | decreasing factor:     %d\n', sensdecfactor0);
    end
    fprintf(' * Mean image     | optimise:              %d\n', meanoptim);
    if meanoptim
        fprintf(' * Mean image     | factor:                %g\n', meanfactorfinal);
        fprintf(' * Mean image     | mean value:            %g\n', meanval);
        fprintf(' * Mean image     | decreasing factor:     %d\n', meandecfactor0);
        fprintf(' * Mean image     | Iter:                  %g\n', meaniter);
        fprintf(' * Mean image     | Iter (initial):        %g\n', meaniter0);
    end
    fprintf(' * Covariance     | optimise:              %d\n', precoptim);
    fprintf(' * Convergence    | tolerance:             %g\n', tol);
    fprintf(' * Convergence    | min iterations:        %d\n', itermin);
    fprintf(' * Convergence    | max iterations:        %d\n', itermax);
    fprintf(' * Convergence    | decreasing iterations: %d\n', iterdec);
    fprintf(' * Input          | voxel size:            [%g %g %g]\n', vs(1), vs(2), vs(3));
    fprintf(' * Input          | acceleration:          %g\n', 1./propmask);
    fprintf(' * Input          | mean:                  %d\n', ~isempty(p.Results.MeanImage));
    fprintf(' * Input          | sensitivities:         %d\n', ~isempty(p.Results.SensMaps));
end

% =========================================================================
%
%                           INITIAL ESTIMATES
%
% =========================================================================

if verbose > -1
    if verbose > 0
        fprintf([repmat('-', [1 78]) '\n']);
    end
    fprintf('Initialisation\n');
end

% -------------------------------------------------------------------------
% Estimate noise precision using Rician fit
% -------------------------------------------------------------------------
% Note that this is the noise precision in k-space, and is therefore 
% independent of the lattice size.
if isempty(prec) || any(any(isnan(prec)))
    if numel(mask) > 1
        % Extract calibration region first
        ac  = utils.ifft(utils.acsub(coils, mask), [1 2 3]);
        Nac = size(ac,1)*size(ac,2)*size(ac,3);
    else
        ac   = coils;
        Nac  = Nvox;
    end
    prec = b1m.init.noise(ac);
    prec = prec / Nac;
    clear ac Nac
end
logdetnoise = utils.logdetPD(prec);


% -------------------------------------------------------------------------
% Compute mask of the background
% -------------------------------------------------------------------------
if isscalar(maskbg)
    if maskbg
        if numel(mask) > 1
            % Extract calibration region first
            ac = utils.ifft(utils.acsub(coils, mask), [1 2 3]);
        else
            ac = coils;
        end
        ac     = sqrt(sum(abs(ac).^2,4));
        dimac  = size(ac);
        ac     = reshape(ac, [], size(ac,5));
        maskbg = utils.gmm.fit(ac, 2);
        maskbg = reshape(maskbg(:,1), dimac);
        spm_diffeo('boundary', 1);
        maskbg = spm_diffeo('resize', single(maskbg), [Nx Ny Nz]);
        maskbg = maskbg > 0.5;
        maskbg = logical(spm_erode(double(maskbg)));
        clear ac dm
    else
        maskbg = [];
    end
end

% -------------------------------------------------------------------------
% Initial estimate of the mean
% -------------------------------------------------------------------------
% One Gauss-Newton update of the mean image, with initial (identity)
% sensitivities
if verbose > 1
    b1m.plot.mean(meanim, prec, ll, vs);
end
if meanoptim && meaninit
    if isnan(lls)
        if sensinitmag && sensinitphase
            lls = zeros(Nc,2);
        else
            % > Initial log-likelihood (sensitivity prior term)
            lls  = b1m.ll.sensitivity(sens, 'RegFactor', sensfactor, 'VoxelSize', vs);
        end
    end
    for i=1:meaniter0
        if verbose > 0, fprintf('> Update mean: '); end
        [meanim,llm,llr,ok,ls] = b1m.update.mean(...
            coils, sens, meanim,             ...
            'RegFactor',     meanfactor,     ...
            'llPrior',       llr,            ...
            'Precision',     prec,           ...
            'VoxelSize',     vs,             ...
            'SensLog',       senslog,        ...
            'BackgroundMask',maskbg,         ...
            'SamplingMask',  mask);
        llm = llm + propmask*Nvox*logdetnoise;
        if verbose > 0
            if ok, fprintf(' :D (%d)\n', ls);
            else,  fprintf(' :(\n');
            end
        end
        ll = [ll (sum(llm)+sum(lls(:))+llr)];
        if verbose > 1
            b1m.plot.mean(meanim, prec, ll, vs);
        end
    end
end

% -------------------------------------------------------------------------
% Initial estimate of the sensitivity
% -------------------------------------------------------------------------
% A few alternated updates of the mean image (Gauss-Newton) and of the 
% (flat) sensitivities. This is to roughly center the mean image (in terms 
% of global magnitude and phase) with respect to the individual coils.
%
% For sensitivity magnitude, the geometric mean, across voxels, of the 
% ratio between the mean image magnitude and each coil image magnitude is
% used.
% 
% For sensitivity phase, the Von Mises mean, across voxels, of the
% difference between the mean image phase and each coil image phase is
% used.
% See: * Bishop's PRML - chapter 2.3.8
%      * https://en.wikipedia.org/wiki/Von_Mises_distribution
if verbose > 1
    b1m.plot.mean(meanim, prec, ll, vs);
end
if sensoptim && sensinit % (sensinitmag || sensinitphase)
    lls = zeros(Nc,2);
    for i=1:3
        if verbose > 0
            fprintf('> Update sensitivity\n');
        end
%         if sensinitmag
%             sens = b1m.init.magnitude(meanim, coils, sens, senslog);
%         end
%         if sensinitphase
%             sens = b1m.init.phase(meanim, coils, sens, senslog);
%         end
        if sensinit
            sens = b1m.init.sensitivity(meanim, coils, sens);
            if senslog
                sens = log(sens);
            end
        end
        
        if sensinit && meanoptim && meaninit
            [sens,meanim] = b1m.centre(sens,meanim);
        end
        
        if meanoptim && meaninit
            if verbose > 0, fprintf('> Update mean: '); end
            [meanim,llm,llr,ok,ls] = b1m.update.mean(...
                coils, sens, meanim,              ...
                'RegFactor',     meanfactor,     ...
                'llPrior',       llr,            ...
                'Precision',     prec,           ...
                'VoxelSize',     vs,             ...
                'SensLog',       senslog,        ...
                'BackgroundMask',maskbg,         ...
                'SamplingMask',  mask);
            llm = llm + propmask*Nvox*logdetnoise;
            if verbose > 0
                if ok, fprintf(' :D (%d)\n', ls);
                else,  fprintf(' :(\n');
                end
            end
            ll = [ll (sum(llm)+sum(lls(:))+llr)];
            if verbose > 1
                b1m.plot.mean(meanim, prec, ll, vs);
            end
        end
end



% -------------------------------------------------------------------------
% Log-likelihood
% -------------------------------------------------------------------------

if isnan(lls)
    % > Initial log-likelihood (sensitivity prior term)
    lls  = b1m.ll.sensitivity(sens, 'RegFactor', sensfactor, 'VoxelSize', vs);
end
if isnan(llm)
    % > Initial log-likelihood (cond term)
    llm = b1m.ll.conditional(coils,sens,meanim,prec,mask,senslog) ...
        + propmask*Nvox*logdetnoise;
end
if isnan(llr)
    % > Initial log-likelihood (mean prior term)
    llr = b1m.ll.mean(meanim, 'RegFactor', meanfactor, 'VoxelSize', vs);
end
ll = [ll (sum(llm)+sum(lls(:))+llr)];

if verbose > 1
    b1m.plot.mean(meanim, prec, ll, vs);
end

% -------------------------------------------------------------------------
% Loop
% -------------------------------------------------------------------------
for it=1:itermax
    
    % ---------------------------------------------------------------------
    % Previous log-likelihood
    % ---------------------------------------------------------------------
    ll0 = ll(end);
    
    % ---------------------------------------------------------------------
    % Update regularisation
    % ---------------------------------------------------------------------
    if numel(sensdecfactor) >= it
        sensfactor = sensfactorfinal .* sensdecfactor(it);
        if it > 1
            lls = lls * sensdecfactor(it) / sensdecfactor(it-1);
        end
    end
    if numel(meandecfactor) >= it
        meanfactor = meanfactorfinal .* meandecfactor(it);
        if it > 1
            llr = llr * meandecfactor(it) / meandecfactor(it-1);
        end
    end
    
    if verbose > -1
        if verbose > 0
            fprintf([repmat('-', [1 78]) '\n']);
        end
        fprintf('Iteration %d\n', it);
        if sensdecfactor0
            if sensoptim
                fprintf('* sens: %7.1e\n', mean(sensfactor(:)));
            end
        end
        if meandecfactor0
            if meanoptim
                fprintf('* mean: %7.1e\n', meanfactor);
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % Update mean/covariance (closed-form)
    % ---------------------------------------------------------------------
%     if optim_prec
%         if verbose > 0
%             fprintf('> Update covariance\n');
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
    if sensoptim
        if verbose > 0, fprintf('> Update Sensitivity:'); end
        if isdiag(prec)
            % If the precision matrix is diagonal, computation of the
            % log-likelihood and derivatives in implemented in a slightly
            % faster way (each sensitivity field only depends on one coil).
            %
            % * This means that we can distribute the processing of each coil!
            %
            % * This also means that we need to sum the log-likelihood of each 
            %   coil to get the complete model log-likelihood.
            lls0 = lls;
            lls  = zeros(Nc,2);
            llm  = zeros(Nc,1);
            Ad   = diag(prec);

            if verbose < 3
                % We can parallelise :D
                parfor(n=1:Nc, double(Nw))
                    [sens(:,:,:,n),llm(n),lls(n,:),ok,ls] = b1m.update.sensitivity(...
                        meanim, coils(:,:,:,n), sens(:,:,:,n),  ...
                        'Precision',     Ad(n),                 ...
                        'RegFactor',     sensfactor(n,:),       ...
                        'VoxelSize',     vs,                    ...
                        'Log',           senslog,               ...
                        'LLPrior',       lls0(n,:),             ...
                        'SamplingMask',  mask);
                    if verbose > 0
                        if ok, fprintf(' [%2d :D (%d)]', n, ls);
                        else,  fprintf(' [%2d :(    ]', n);
                        end
                    end
                end
            else
                % We cannot parallelise :(
                for n=1:Nc
                    if verbose > 0
                        if verbose > 2
                            b1m.plot.fit(n, coils, sens, meanim, mask, senslog, vs);
                        end
                        if mod(n,4) == 1, fprintf('\n  | '); end
                        fprintf('%2d', n);
                    end
                    [sens(:,:,:,n),llm(n),lls(n,:),ok,ls] = b1m.update.sensitivity(...
                        meanim, coils(:,:,:,n), sens(:,:,:,n),  ...
                        'Precision',     Ad(n),                 ...
                        'RegFactor',     sensfactor(n,:),       ...
                        'VoxelSize',     vs,                    ...
                        'Log',           senslog,               ...
                        'LLPrior',       lls0(n,:),             ...
                        'SamplingMask',  mask);
                    if verbose > 0
                        if ok, fprintf(' :D (%d) | ', ls);
                        else,  fprintf(' :(     | ');
                        end
                        if verbose > 2
                            b1m.plot.fit(n, coils, sens, meanim, mask, senslog, vs);
                        end
                    end
                end
            end
            llm = sum(llm);
        else
            % We do not parallelise, but still need to sum individual 
            % coil-wise prior log-likelihoods
            lls0 = lls;
            lls  = zeros(Nc,2);
            llm  = 0;
            for n=1:Nc

                if verbose > 0
                    if verbose > 2
                        b1m.plot.fit(n, coils, sens, meanim, mask, senslog, vs);
                    end
                    if mod(n,4) == 1, fprintf('\n  | '); end
                    fprintf('%2d', n);
                end
                [sens,llm,lls,ok,ls] = b1m.update.sensitivity(...
                    meanim, coils, sens,                ...
                    'Index',         n,                 ...
                    'Precision',     prec,              ...
                    'RegFactor',     sensfactor,        ...
                    'VoxelSize',     vs,                ...
                    'Log',           senslog,           ...
                    'LLPrior',       lls0,              ...
                    'SamplingMask',  mask);
                if verbose > 0
                    if ok(n), fprintf(' :D (%d) | ', ls);
                    else,     fprintf(' :(     | ');
                    end
                    if verbose > 2
                        b1m.plot.fit(n, coils, sens, meanim, mask, senslog, vs);
                    end
                end

            end
        end
        if verbose > 0, fprintf('\n'); end
        % Add normalisation term
        llm = llm + propmask*Nvox*logdetnoise;
        ll  = [ll (sum(llm)+sum(lls(:))+llr)];
        if verbose > 1
            b1m.plot.mean(meanim, prec, ll, vs);
        end
    end
    
    % ---------------------------------------------------------------------
    % Centre
    % ---------------------------------------------------------------------
    if senslog
        [sens,meanim] = b1m.centre(sens, meanim, sensfactor);
        if verbose > 1
            b1m.plot.mean(meanim, prec, ll, vs);
        end
    end
    
    % ---------------------------------------------------------------------
    % Update mean (Gauss-Newton)
    % ---------------------------------------------------------------------
    if meanoptim
        if verbose > 0, fprintf('> Update mean:\n '); end
        for i=1:1:meaniter
            [meanim,llm,llr,ok,ls] = b1m.update.mean(...
                coils, sens, meanim,             ...
                'RegFactor',     meanfactor,     ...
                'llPrior',       llr,            ...
                'Precision',     prec,           ...
                'VoxelSize',     vs,             ...
                'SensLog',       senslog,        ...
                'BackgroundMask',maskbg,         ...
                'SamplingMask',  mask);
            llm = llm + propmask*Nvox*logdetnoise;
            if verbose > 0
                if ok, fprintf(' :D (%d) |', ls);
                else,  fprintf(' :(     |');
                end
            end
            ll = [ll (sum(llm)+sum(lls(:))+llr)];
            if verbose > 1
                b1m.plot.mean(meanim, prec, ll, vs);
            end
       end
       if verbose > 0, fprintf('\n'); end
    end
    
%     % ---------------------------------------------------------------------
%     % Update log-likelihood
%     % ---------------------------------------------------------------------
%     ll = [ll (sum(llm)+sum(lls(:))+llr)];
%     if verbose > 1
%         b1m.plot.mean(meanim, prec, ll, vs);
%     end
    
    % ---------------------------------------------------------------------
    % Check gain
    % ---------------------------------------------------------------------
    if it > 1
        llmin = min(ll, [], 'omitnan');
        llmax = max(ll, [], 'omitnan');
        gain  = (ll(end) - ll0)/(llmax-llmin);
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

if isreal(coils)
    sens = real(sens);
end

% -------------------------------------------------------------------------
% Time execution
% -------------------------------------------------------------------------
if verbose > -1
    stop = toc(start);
    fprintf('Processing finished: in %s\n', utils.sec2ydhms(stop));
end

end % < function multicoil_infer