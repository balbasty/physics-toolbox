function [dy,info] = relax(H, g, w, mode, prec, vs, opt, dy)
% Relaxation solver for L1 spatial regularisation.
%
% FORMAT d = mpm.l1.solver.cg(H, g, w, mode, prec, vs, [opt])
% H    - {nx ny nz nf(nf+1)/2} - Field of (sparse) Hessian matrices
% g    - {nx ny nz nf}         - Field of gradients
% w    - {nx ny nz}            - Weights (reweighted Gaussian approximation)
% mode - {nf}                  - Regularisation mode for each channel (0|1|2)
% prec - {nf}                  - Regularisation precision for each channel
% vs   - {3}                   - Voxel size
% opt  - {struct}              - Structure of options with fields:
%      . nbiter                - Number of CG iterations [10]
%      . tolerance             - Early stopping tolerance [1E-3]
%      . verbose               - Print progress [true]
% d    - {nx ny nz nf}         - Step: d = H\g

    if nargin < 8, dy   = []; end
    if nargin < 7, opt  = struct; end
    if ~isfield(opt, 'fmginit'),   opt.fmginit   = true;  end

    % Initial guess using a majoriser of the true Hessian
    wbnd = double(max(w(:)));
    prec_bound = prec;
    prec_bound(mode==1) = wbnd * prec_bound(mode==1);
    spm_field('boundary', 1);
    start = tic;
    if isempty(dy)
        if opt.fmginit
            dy = spm_field(H, g, [vs 0 1 0 2 2], prec_bound(:)');
        else
            dy = zeros(size(g), 'single');
        end
    end
    time_init = toc(start);
    
    if all(w(:)==1)
        % No need to relax
        return
    end
    
    % Smoothing term: use diagonal of membrane kernel
    spm_diffeo('boundary', 1);
    scl = spm_diffeo('kernel', [3 3 3], [vs 0 1 0 0 0]);
    scl = double(abs(scl(1,1,1)));
    scl = scl*wbnd;
    
    % Define forward and inverse operators
    function Ax = A(x,ind)
    % Forward operator: A = E + F
        spm_field('boundary', 1);
        if nargin < 2
            Ax = zeros(size(x), 'like', x);
            Ax(:,:,:,mode==1) = spm_field('vel2mom1', single(x(:,:,:,mode==1)), single(w), [vs   1  ], prec(mode==1));
            Ax(:,:,:,mode==2) = spm_field('vel2mom',  single(x(:,:,:,mode==2)),            [vs 0 1 0], prec(mode==2));
            Ax = Ax + spm_field('Atimesp', H, x);
        else
            Ax = zeros(size(x), 'like', x);
            Ax(:,:,:,mode==1) = spm_field('vel2mom1', single(x(:,:,:,mode==1)), single(w), [vs   1  ], prec(mode==1));
            Ax(:,:,:,mode==2) = spm_field('vel2mom',  single(x(:,:,:,mode==2)),            [vs 0 1 0], prec(mode==2));
            Ax   = reshape(Ax, [], size(Ax,4));
            Ax   = Ax(ind,:);
            Ax   = reshape(Ax, [], 1, 1, size(Ax,2));
            Hsub = reshape(H, [], size(H,4));
            Hsub = Hsub(ind,:);
            Hsub = reshape(Hsub, [], 1, 1, size(Hsub,2));
            xsub = reshape(x, [], size(x,4));
            xsub = xsub(ind,:);
            xsub = reshape(xsub, [], 1, 1, size(xsub,2));
            Ax   = Ax + spm_field('Atimesp', Hsub, xsub);
            Ax   = reshape(Ax, [], size(Ax,4));
        end
    end
    function x = iE(x,ind)
    % Forward operator: inv(E)
        spm_field('boundary', 1);
        if nargin < 2
            x = spm_field(H, x, [1 1 1 scl 0 0], prec);
        else
            Hsub = reshape(H, [], size(H,4));
            Hsub = Hsub(ind,:);
            Hsub = reshape(Hsub, [], 1, 1, size(Hsub,2));
            xsub = reshape(x, [], size(x,4));
            xsub = xsub(ind,:);
            xsub = reshape(xsub, [], 1, 1, size(xsub,2));
            x    = spm_field(Hsub, xsub, [1 1 1 scl 0 0], prec);
            x    = reshape(x, [], size(x,4));
        end
        
    end
    
    % Relax
    [dy,info] = optim.relax3(@A, g, @iE, dy, opt);
    info.time_init = time_init;
    
end

