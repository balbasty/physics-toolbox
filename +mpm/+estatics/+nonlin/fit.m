function out = fit(in,opt)
% Nonlinear ESTATICS fit.
%
%   The ESTATICS model applies to multi-echo spoiled gradient-echo 
%   acquisitions. It assumes that all echo series possess the same R2* 
%   decay, but different initial signals (at TE=0).
%   This function fits the ESTATICS model using [regularised] non-linear
%   weighted least-squares.
%
% FORMAT out = mpm.estatics.nonlin.fit(in,opt)
%
% in  - Structure of input data generated by `mpm.io.input`
%
% opt - Structure of parameters with (optional) fields:
%   nbscales     [5]       Number of scales (or zoom levels) 
%   nbiter       [3]       Number of iterations per scale
%   tolerance    [1E-4]    Gain threshold for early stopping (per scale)
%   out.folder   ['.']     Output folder 
%   out.fname    ['.nii']  Basis for output filenames 
%   out.mem      ['map']   Memory map or load output volumes (map|load)
%   reg.mode     [0 1]     Regularisation type (0=none|1=l1/tv|2=l2)
%   reg.prec.<c> [10 5E3]  Regularisation factor
%           .R2s [10 5E0]
%   reg.mean.<c> [NaN]     Prior mean
%   fov          [0]       Index of input FOV used for model space (0=bbox)
%   vs           [NaN]     Voxel size of model space
%   coreg        [true]    Initial co-registration
%   init         ['mean']  Init mode: 'mean'/'logfit'/'minilogfit'
%   subsample    [Inf]     Subsampling distance in mm
%   threads      [-1]      Number of threads used by SPM
%   solver.type      ['relax'] Solver used for L1 regularisation (relax|cg)
%   solver.nbiter    [10]      Number of iterations of the linear solver
%   solver.tolerance [0]       Solver gain threshold for early stopping
%   solver.verbose   [1]       Solver verbosity 
%
% Usually: `NaN` means 'informed default'
%          `Inf` means 'as precise as possible'
% reg.mode and reg.prec take two values: (1) absolute valuefPD,fT1,fMTs
%                                        (2) spatial gradients


% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
if nargin < 2, opt = struct; end
opt = mpm.estatics.nonlin.opt(opt);

% -------------------------------------------------------------------------
% Multithread SPM
% -------------------------------------------------------------------------
nthreads0 = utils.threads.get();
utils.threads.set(opt.threads);

% -------------------------------------------------------------------------
% Index available contrasts (only keep FLASH volumes)
% -------------------------------------------------------------------------
issge       = cellfun(@(x) strcmpi(x.seq, 'SGE'), in);
in          = in(issge);
types       = cellfun(@(x) x.type, in, 'UniformOutput', false);
contrasts   = [unique(types) {'R2s'}];
nbcontrasts = numel(contrasts);

% -------------------------------------------------------------------------
% Default regularisation parameters
% -------------------------------------------------------------------------
% --- set default parameters
for ct=contrasts
    ct = ct{1};
    if~isfield(opt.reg.prec, ct)
        opt.reg.prec.(ct) = opt.reg.prec.default;
    end
    if~isfield(opt.reg.mean, ct)
        opt.reg.mean.(ct) = opt.reg.mean.default;
    end
    if~isfield(opt.reg.mode, ct)
        opt.reg.mode.(ct) = opt.reg.mode.default;
    end
end
% --- convert named parameters to vectors
reg0 = opt.reg;
opt.reg.mode = cell2mat(cellfun(@(x) reg0.mode.(x), contrasts(:), 'UniformOutput', false));
opt.reg.prec = cell2mat(cellfun(@(x) reg0.prec.(x), contrasts(:), 'UniformOutput', false));
opt.reg.mean = cellfun(@(x) reg0.mean.(x), contrasts(:));
% --- estimate mean parameter value + correct regularisation
opt.verbose > 0 && fprintf('Guess mean parameter value to correct regularisation\n');
mu = mpm.estatics.loglin.fit.mini(in);
mu = cellfun(@(x) mu.(x).dat, contrasts(:));
opt.reg.mean(isnan(opt.reg.mean)) = mu(isnan(opt.reg.mean));
opt.reg.prec(:,2) = opt.reg.prec(:,2) ./ abs(mu(:));
for n=1:numel(mu)
    opt.verbose > 0 && fprintf('* Param %d | mean: %7.3g | prec: %7.3g\n', n, mu(n), opt.reg.prec(n,2));
end
% TODO: ^ is this correct for L2 regularisation?

% -------------------------------------------------------------------------
% Co-registration
% -------------------------------------------------------------------------
opt.verbose > 0 && fprintf('Coregister volumes\n');
if opt.coreg
    in = mpm.coregister(in);
end

% -------------------------------------------------------------------------
% Prepare output
% -------------------------------------------------------------------------
opt.verbose > 0 && fprintf('Prepare output data structure\n');
% --- Field of view
dim = zeros(3,0);
mat = zeros(4, 4, 0);
for v=1:numel(in)
    dim(:,end+1)   = in{v}.dim(:)';
    mat(:,:,end+1) = in{v}.mat(:,:);
end
[dim,mat] = utils.fov(dim, mat, opt.vs, opt.fov);
% --- Multiscale settings
scales = mpm.compute_scales(opt.nbscales, dim, mat, opt.subsample);
% --- Prefix/Value
prefix = contrasts;
value  = zeros(1,nbcontrasts);
dim    = repmat([scales(end).dim(:); 1], [1 nbcontrasts]);
mat    = repmat(scales(end).mat, [1 1 nbcontrasts]);
if any(opt.reg.mode(:,2)==1)
    prefix = [prefix {'W'}];
    value  = [value 1];
    dim    = [dim [scales(end).dim(:); 1]];
    if ischar(opt.reg.uncertainty) && strcmpi(opt.reg.uncertainty, 'bayes')
        prefix = [prefix {'U'}];
        value  = [value 1E-3];
        dim    = [dim [scales(end).dim(:); nbcontrasts]];
    end
end
% --- Create
out = mpm.io.output(prefix, dim, mat, value, 'single', opt.out);
% --- Save TR/FA
for v=1:numel(in)
    out.(in{v}.type).extras.TR = in{v}.TR;
    out.(in{v}.type).extras.FA = in{v}.FA;
end

% -------------------------------------------------------------------------
% Initialise with loglinear fit
% -------------------------------------------------------------------------
switch lower(opt.init)
    case 'logfit'
        opt.verbose > 0 && fprintf('Initial loglinear fit\n');
        out = mpm.estatics.loglinfit.core(in,out,opt.subsample,opt.verbose);
    case 'mean'
        for k=1:numel(contrasts)
            out.(contrasts{k}).dat(:) = opt.reg.mean(k);
        end
    otherwise
        error('Initialisation mode ''%s'' not implemented.', init);
end
opt.verbose > 1 && mpm.estatics.plot.progress(out,[]);

% -------------------------------------------------------------------------
% Nonlinear fit
% -------------------------------------------------------------------------
opt.verbose > 0 && fprintf('Full nonlinear fit\n');
ll = [];
opt.reg.prec0 = opt.reg.prec;
for s=numel(scales):-1:1
    
    % ---------------------------------------------------------------------
    % Get scale specific parameters
    % ---------------------------------------------------------------------
    opt.verbose > 0 && fprintf('Scale %i\n', s);
    dim = scales(s).dim;
    % vs  = scales(s).vs;
    vs  = sqrt(sum(scales(s).mat(1:3,1:3).^2));
    sub = scales(s).sub;
    ff  = scales(s).ff;
    opt.reg.prec = ff * opt.reg.prec0 * prod(vs);
    
    % ---------------------------------------------------------------------
    % A bit of ad-hoc fudging for L2 regularisation
    % ---------------------------------------------------------------------
    opt.reg.prec(opt.reg.mode(:,2)==2,2) = opt.reg.prec(opt.reg.mode(:,2)==2,2) * power(10,s-1);
    
    it0 = numel(ll)+1;
    for it=1:(opt.nbiter+s)
    
        opt.verbose > 0 && fprintf('Iteration %i\n', it);

        % -----------------------------------------------------------------
        % Update maps
        % -----------------------------------------------------------------
        [ind,K] = utils.symIndices(nbcontrasts);   % Sparse indices for symmetric matrices
        g   = zeros([dim nbcontrasts], 'single');  % Gradient
        H   = zeros([dim K], 'single');            % Hessian (symmetric 4x4 matrix)
        llx = 0;                                   % Log-likelihood: data term
        lly = 0;                                   % Log-likelihood: prior term

        % -----------------------------------------------------------------
        % Loop over volumes
        opt.verbose > 0 && fprintf('Gradient: data ');
        for v=1:numel(in)
            opt.verbose > 0 && fprintf('%d',v);

            % - Compute gradient
            [llx1,g1,H1] = mpm.estatics.nonlin.gradient(in{v}, out, sub);
            llx = llx + llx1;

            % - Get index
            k = find(strcmpi(in{v}.type, contrasts));   % Intercept
            l = find(strcmpi('R2s', contrasts));        % Decay
            
            % - Add to full gradient
            gind = [k l];
            g(:,:,:,gind) = g(:,:,:,gind) + g1;
            clear g1
            Hind = [ind(k,k) ind(l,l) ind(k,l)];
            H(:,:,:,Hind) = H(:,:,:,Hind) + H1;
            clear H1
            
            opt.verbose > 0 && fprintf(' ');
        end
        opt.verbose > 0 && fprintf('\n');

        % -----------------------------------------------------------------
        % Gradient: Absolute
        if any(opt.reg.mode(:,1) > 0)
            opt.verbose > 0 && fprintf('Gradient: absolute ');
            for k=1:nbcontrasts
                if opt.reg.mode(k,1) == 2
                    fprintf('.');
                    y = single(out.(contrasts{k}).dat()) - opt.reg.mean(k);
                    g(:,:,:,k) = g(:,:,:,k) + opt.reg.prec(k,1) * y;
                    H(:,:,:,k) = H(:,:,:,k) + opt.reg.prec(k,1);
                    lly = lly - 0.5 * sum(y(:).^2, 'double');
                    clear y
                end
            end
            opt.verbose > 0 && fprintf('\n');
        end

        % -----------------------------------------------------------------
        % Gradient: Membrane
        if any(opt.reg.mode(:,2) > 0)
            opt.verbose > 0 && fprintf('Gradient: membrane ');
            w    = 0;
            wnew = 0;
            if any(opt.reg.mode(:,2) == 1), w = 1./single(out.W.dat()); end
            for k=1:nbcontrasts
                switch opt.reg.mode(k,2)
                    case 1
                        opt.verbose > 0 && fprintf('.');
                        Ly = single(out.(contrasts{k}).dat());
                        [Ly,Dy] = mpm.l1.vel2mom(Ly, opt.reg.prec(k,2), vs, w);
                        Dy = sum(sum(Dy.^2,5),4);
                        wnew = wnew + Dy;
                        g(:,:,:,k) = g(:,:,:,k) + Ly;
                        clear Ly
                    case 2
                        opt.verbose > 0 && fprintf('.');
                        y  = single(out.(contrasts{k}).dat());
                        Ly = mpm.l2.vel2mom(y, opt.reg.prec(k,2), vs);
                        g(:,:,:,k) = g(:,:,:,k) + Ly;
                        lly = lly - 0.5*sum(y(:).*Ly(:), 'double');
                        clear y Ly
                end
            end
            opt.verbose > 0 && fprintf('\n');
            if any(opt.reg.mode(:,2) == 1)
                if ischar(opt.reg.uncertainty) && strcmpi(opt.reg.uncertainty, 'bayes')
                    wnew = sqrt(wnew + sum(single(out.U.dat()), 4));
                else
                    wnew = sqrt(wnew + opt.reg.uncertainty);
                end
                lly = lly - sum(wnew(:), 'double');
            end
        end

        % -----------------------------------------------------------------
        % Update MTV weights
        if any(opt.reg.mode(:,2)==1)
            opt.verbose > 0 && fprintf('Update: MTV weights\n');
            out.W.dat(:,:,:) = wnew; clear Wnew
        end

        % -----------------------------------------------------------------
        % Log-likelihood
        ll = [ll llx+lly];

        % -----------------------------------------------------------------
        % Load diagonal of the Hessian
        H(:,:,:,1:nbcontrasts) = bsxfun(@plus, H(:,:,:,1:nbcontrasts), sum(H(:,:,:,1:nbcontrasts)) * eps('single'));

        % -----------------------------------------------------------------
        % Gauss-Newton
        % - Compute step
        if all(opt.reg.mode(:,2) == 0)
        % No regularisation
            dy = mpm.l0.solve(H,g);
        elseif any(opt.reg.mode(:,2) == 1)
        % Some L1 regularisation
            dy = opt.solver.fun(H, g, w, opt.reg.mode(:,2), opt.reg.prec(:,2), vs, opt.solver);
            % Uncertainty about y
            if ischar(opt.reg.uncertainty) && strcmpi(opt.reg.uncertainty, 'bayes')
                out.U.dat(:,:,:,:) = mpm.l1.uncertainty(H, opt.reg.prec(:,2), vs, w);
            end
        else
        % Some L2 but no L1 regularisation
            dy = mpm.l2.solve(H, g, opt.reg.mode(:,2), opt.reg.prec(:,2), vs);
        end
        clear H g W
        % - Update
        opt.verbose > 0 && fprintf('Update: maps\n');
        for k=1:nbcontrasts
            ct = contrasts{k};
            switch ct
                case 'R2s'
                    out.(ct).dat(:,:,:) = max(0, out.(ct).dat(:,:,:) - dy(:,:,:,k));
                otherwise
                    out.(ct).dat(:,:,:) = out.(ct).dat(:,:,:) - dy(:,:,:,k);
            end
        end
        clear dy

        % -----------------------------------------------------------------
        % Gain
        if numel(ll) > it0
            gain = abs((ll(end) - ll(end-1)) / (ll(end-1) - min(ll(it0:end))));
        else
            gain = Inf;
        end

        % -----------------------------------------------------------------
        % Plot
        opt.verbose > 1 && mpm.estatics.plot.progress(out,ll);
        opt.verbose > 0 && fprintf('%s\n', repmat('-',[1 80]));
        opt.verbose > 0 && fprintf('ll = %7.3g | llx = %7.3g | lly = %7.3g | gain = %7.3g\n', ll(end), llx, lly, gain);
        opt.verbose > 0 && fprintf('%s\n', repmat('-',[1 80]));
       
        % -----------------------------------------------------------------
        % Out?
        if gain < opt.tolerance
            break
        end

    end
        
    % ---------------------------------------------------------------------
    % Resize maps
    % ---------------------------------------------------------------------
    if s ~= 1
        out = mpm.resize_output(out, scales(s-1).dim);
    end
    
end

out = mpm.estatics.extrapolate(out,opt.out.fname);

utils.threads.set(nthreads0);