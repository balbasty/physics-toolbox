function [out,ll] = core(in,out,opt)
% Loglinear ESTATICS fit [core algorithm]
%
%   The ESTATICS model applies to multi-echo spoiled gradient-echo 
%   acquisitions. It assumes that all echo series possess the same R2* 
%   decay, but different initial signals (at TE=0).
%   This function fits the ESTATICS model using least-squares in the log
%   domain.
%
% This function only implements the core algorithm and expects the input
% and output structures to be already prepared. It is advised to use the
% high-level wraper `mpm.estatics.loglin.fit` instead.
%
% FORMAT out = mpm.estatics.loglin.core(in,out,,[opt])
% in  - Input data structure (see `mpm.io.input`)
% out - Output data structure (see `mpm.io.output`)
% opt - Structure of parameters with [optional] fields:
%   . scaled    - Scale noise precision by squared signal [false]
%   . subsample - Subsampling distance in mm [Inf=no]
%   . opt.verbose   - Verbosity (0=quiet|[1]=print|2=plot)

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
if nargin < 3, opt = struct; end
if ~isfield(opt, 'subsample'), opt.subsample = Inf;   end
if ~isfield(opt, 'verbose'),   opt.verbose   = 0;     end
if ~isfield(opt, 'scaled'),    opt.scaled    = false; end

% -------------------------------------------------------------------------
% Get avaiable contrasts
% -------------------------------------------------------------------------
issge       = cellfun(@(x) strcmpi(x.seq, 'SGE'), in);
in          = in(issge);
types       = cellfun(@(x) x.type, in, 'UniformOutput', false);
contrasts   = [unique(types) {'R2s'}];
nbcontrasts = numel(contrasts);

% -------------------------------------------------------------------------
% Loglinear fit
% -------------------------------------------------------------------------
ll = [];
for it=1:3
    
    opt.verbose > 0 && fprintf('Iteration %i\n', it);
    
    % ---------------------------------------------------------------------
    % Update maps
    % ---------------------------------------------------------------------
    [ind,K] = symIndices(nbcontrasts);              % Sparse indices for symmetric matrices
    g   = zeros([out.dim nbcontrasts], 'single');   % Gradient
    H   = zeros([out.dim K], 'single');             % Hessian (symmetric 4x4 matrix)
    llx = 0;                                        % Log-likelihood: data term
    
    % ---------------------------------------------------------------------
    % Loop over echoes (can be parallelised using parfor)
    opt.verbose > 0 && fprintf('Gradient: data ');
    for v=1:numel(in)
        opt.verbose > 0 && fprintf('%d',v);

        % - Compute gradient
        [llx1,g1,H1] = mpm.estatics.loglin.gradient(in{v}, out, opt);
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
    
    
    % ---------------------------------------------------------------------
    % Log-likelihood
    ll = [ll llx];
    
    % ---------------------------------------------------------------------
    % Load diagonal of the Hessian
    H(:,:,:,1:nbcontrasts) = bsxfun(@plus, H(:,:,:,1:nbcontrasts), sum(H(:,:,:,1:nbcontrasts)) * eps('single'));
    
    % ---------------------------------------------------------------------
    % Gauss-Newton
    opt.verbose > 0 && fprintf('Gauss-Newton update\n');
    dy = mpm.l0.solve(H,g);
    clear H g
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
    
    % ---------------------------------------------------------------------
    % Plot
    opt.verbose > 1 && mpm.estatics.plot.progress(out,ll);
end