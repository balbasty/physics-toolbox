function [mpm,in] = mpm_nonlinfit(fPD,fT1,fMT,fB1t,fB1r,opt)

if nargin < 6 || isempty(opt)
    opt = struct;
end
if nargin < 5 || isempty(fB1r)
    fB1r = {};
end
if nargin < 4 || isempty(fB1t)
    fB1t = {};
end
if nargin < 3 || isempty(fMT)
    fMT = {};
end
if nargin < 2 || isempty(fT1)
    fT1 = {};
end
opt.subsample = [Inf Inf Inf]; % 

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
itermin   = 20;              % Minimum number of iterations
itermax   = 20;             % Maximum number of iterations
threshold = 1E-4;           % Gain threshold for early stopping
outfolder = './output_mpm/';    % Output folder
outfname  = 'test.nii';     % Prefix for output files 
reg_mode  = [2 1];          % Absolute/Membrane regul.: 0=None, 1=L1/MTV, 2=L2
reg_prec  = [1E-1 10];          % Absolute/Membrane precision (NaN = auto)
reg_mean  = [9; 3; 0; -4];          % Mean (log) parameter value(NaN = auto)
vs        = [0.8 0.8 0.8];        % Reconstruction voxel size (Nan = from input)
if isfield(opt, 'itermin'),   itermin   = opt.itermin;   end
if isfield(opt, 'itermax'),   itermax   = opt.itermax;   end
if isfield(opt, 'threshold'), threshold = opt.threshold; end
if isfield(opt, 'outfolder'), outfolder = opt.outfolder; end
if isfield(opt, 'outfname'),  outfname  = opt.outfname;  end
if isfield(opt, 'reg_mode'),  reg_mode  = opt.reg_mode;  end
if isfield(opt, 'reg_prec'),  reg_prec  = opt.reg_prec;  end
if isfield(opt, 'reg_mean'),  reg_mean  = opt.reg_mean;  end
if isfield(opt, 'subsample'), subsample = opt.subsample; end
if isfield(opt, 'vs'),        vs        = opt.vs; end

reg_mode = padutil(reg_mode, [4-size(reg_mode,1) 0], 'replicate', 'post');
reg_prec = padutil(reg_prec, [4-size(reg_prec,1) 0], 'replicate', 'post');
reg_mean = padutil(reg_mean, [4-size(reg_mean,1) 0], 'replicate', 'post');

% reg_prec(:,2) = reg_prec(:,2) * prod(vs)^(1/3);

% -------------------------------------------------------------------------
% Read input
% -------------------------------------------------------------------------
fprintf('Prepare input data structure\n');
in = prepare_input(fPD, fT1, fMT);

% -------------------------------------------------------------------------
% Co-registration
% -------------------------------------------------------------------------
fprintf('Coregister volumes\n');
in = coregister(in);

% -------------------------------------------------------------------------
% Read precomputed bias field
% -------------------------------------------------------------------------
pre = prepare_precomputed(fB1t, fB1r);
% TODO: apply transformation matrices
for v=1:size(in.trf, 3)
    i_B1t = find(pre.B1t.vol == v);
    if isempty(i_B1t)
        i_B1t = find(pre.B1t.vol == Inf);
    end
    pre.B1t.mat(:,:,i_B1t) = in.trf(:,:,v)\pre.B1t.mat0(:,:,i_B1t);
end
for v=1:size(in.trf, 3)
    i_B1r = find(pre.B1r.vol == v);
    if isempty(i_B1r)
        i_B1r = find(pre.B1r.vol == Inf);
    end
    pre.B1r.mat(:,:,i_B1r) = in.trf(:,:,v)\pre.B1r.mat0(:,:,i_B1r);
end

% -------------------------------------------------------------------------
% Prepare output
% -------------------------------------------------------------------------
fprintf('Prepare output data structure\n');
[odim,omat] = compute_fov(in.dim,in.mat,vs);
odim        = odim(:)'; % ensure row vector
prefix      = {'PDw' 'T1w' 'MTw' 'R2'};
value       = [0 0 0 1];
if any(reg_mode(:,2)==1)
    prefix  = [prefix {'W'}];
    value   = [value 1];
end
% out = prepare_output(prefix, odim, omat, outfolder, ['estatics_' outfname], value);

% -------------------------------------------------------------------------
% Initialise with loglinear fit
% -------------------------------------------------------------------------
% fprintf('Initial loglinear fit + rational approximation\n');
% out = estatics_loglinfit(in,out,opt);
% mpm = estatics_to_mpm(out,in,pre,true);
% if isfield(out, 'W'), mpm.W = out.W; end

fprintf('Initialise maps\n');
mpm = prepare_output({'A' 'R2' 'R1' 'MT' 'W'}, odim, omat, outfolder, 'mpm.nii', [reg_mean(:)' 1]);
mpm.dat = cat(4, mpm.A.dat, mpm.R2.dat, mpm.R1.dat, mpm.MT.dat);
mpm.A.idx  = 1;
mpm.R2.idx = 2;
mpm.R1.idx = 3;
mpm.MT.idx = 4;

mpm_plot_progress(mpm,[]);

% -------------------------------------------------------------------------
% Nonlinear fit
% -------------------------------------------------------------------------
fprintf('Full nonlinear fit\n');
ll = [];
decfactor = logspace(4,0,8);
reg_prec0 = reg_prec;
for it=1:itermax
    
    fprintf('Iteration %i\n', it);
    
    reg_prec = reg_prec0;
    if ~isempty(decfactor)
        df = decfactor(1);
        decfactor = decfactor(2:end);
        fprintf('Factor: %g\n', df);
        reg_prec(:,1) = reg_prec(:,1) * df;
        df = sqrt(df);
        reg_prec(:,2) = (reg_prec(:,2) * df) .* (df * (reg_mode(:,2) == 1));
    end
    
    % ---------------------------------------------------------------------
    % Update maps
    % ---------------------------------------------------------------------
    K = 4;
    [~,KK] = symIndices(K);            % Sparse indices for symmetric matrices
    g   = zeros([odim K],  'single');  % Gradient
    H   = zeros([odim KK], 'single');  % Hessian (symmetric 4x4 matrix)
    llx = 0;                           % Log-likelihood: data term
    lly = 0;                           % Log-likelihood: prior term
    
    % ---------------------------------------------------------------------
    % Loop over echoes (can be parallelised using parfor)
    fprintf('Gradient: data ');
    parfor n=1:numel(in.dat)
        fprintf('.');
        
        [llx1,g1,H1] = gradient_mpm_nonlin(n, in, mpm, pre, opt);
        
        llx = llx + llx1;
        g   = g   + g1;
        H   = H   + H1;
    end
    fprintf('\n');
    
    % ---------------------------------------------------------------------
    % Gradient: Absolute
    if any(reg_mode(:,1) > 0)
        fprintf('Gradient absolute ');
        parfor k=1:4
            if reg_mode(k,1) == 2
                fprintf('.');
                y = single(mpm.dat(:,:,:,k)) - reg_mean(k);
                g(:,:,:,k) = g(:,:,:,k) + reg_prec(k,1) * y;
                H(:,:,:,k) = H(:,:,:,k) + reg_prec(k,1);
                lly = lly - 0.5 * reg_prec(k,1) * sum(y(:).^2);
                y = [];
            end
        end
        fprintf('\n');
    end
    
    % ---------------------------------------------------------------------
    % Gradient: Membrane
    if any(reg_mode(:,2) > 0)
        fprintf('Gradient membrane ');
        W    = 0;
        Wnew = 0;
        if any(reg_mode(:,2) == 1)
            W    = mpm.W.dat();
        end
        parfor k=1:4
            if reg_mode(k,2) == 1
                fprintf('.');
                y  = single(mpm.dat(:,:,:,k));
                [Ly,Dy] = vel2mom_l1(y, reg_prec(k,2), vs, W);
                y  = [];
                Dy = sum(sum(Dy.^2,5),4);
                g(:,:,:,k) = g(:,:,:,k) + Ly;
                Ly = [];
                Wnew = Wnew + Dy;
                Dy = [];
            elseif reg_mode(k,2) == 2
                fprintf('.');
                y  = single(mpm.dat(:,:,:,k));
                Ly = vel2mom_l2(y, reg_prec(k,2), vs);
                g(:,:,:,k) = g(:,:,:,k) + Ly;
                lly = lly - 0.5*y(:)'*Ly(:);
                y = [];
                Dy = [];
            end
        end
        fprintf('\n');
        if any(reg_mode(:,2) == 1)
            Wnew = sqrt(abs(Wnew) + 1e-5);
            lly = lly - sum(Wnew(:));
        end
    end
    
    % ---------------------------------------------------------------------
    % Update MTV weights
    if any(reg_mode(:,2)==1) && it >= 12
        fprintf('Update MTV weights\n');
        mpm.W.dat(:,:,:) = Wnew; clear Wnew
    end
    
    % ---------------------------------------------------------------------
    % Log-likelihood
    ll = [ll llx+lly];
    
    % ---------------------------------------------------------------------
    % Load diagonal of the Hessian
    H(:,:,:,1:4) = bsxfun(@plus, H(:,:,:,1:4), max(H(:,:,:,1:4)) * 1e-4);
    
    % ---------------------------------------------------------------------
    % Gauss-Newton
    fprintf('Gauss-Newton update\n');
    if all(reg_mode(:,2) == 0)
        dy = solve_l0(H,g);
    elseif any(reg_mode(:,2) == 1)
        dy = solve_l1(H, g, reg_mode(:,2), reg_prec(:,2), vs, W);
    else % Some L2 but no L1 regularisation
        dy = solve_l2(H, g, reg_mode(:,2), reg_prec(:,2), vs);
    end
    clear H g
    
    mpm.dat(:,:,:,:) = mpm.dat() - 0.5 * dy;
    clear dy
    
    % ---------------------------------------------------------------------
    % Gain
    if numel(ll) > 1
        gain = abs((ll(end) - ll(end-1)) / (ll(end-1) - min(ll)));
    else
        gain = Inf;
    end
        
    % ---------------------------------------------------------------------
    % Plot
    mpm_plot_progress(mpm,ll);
    fprintf('ll = %7.3g | llx = %7.3g | lly = %7.3g | gain = %7.3g\n', ll(end), llx, lly, gain);
    
    % ---------------------------------------------------------------------
    % Out?
    if it >= itermax || (it >= itermin && gain < threshold)
        break
    end
    
end

% -------------------------------------------------------------------------
% Log to proper unit
for k=1:size(mpm.dat,4)
    mpm.dat(:,:,:,k) = exp(mpm.dat(:,:,:,k));
end
if size(mpm.dat,4) >= 4
    mpm.dat(:,:,:,k) = 100 * mpm.dat(:,:,:,k);
end