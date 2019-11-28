function [llx,g,H,res] = gradient_mpm_nonlin(n, in, out, pre, opt)
% FORMAT [llx,g,H,r] = gradient_mpm_nonlin(n, in, out, [pre], [opt])
%
% INPUT
% -----
% n   - Index of input volume to process
% in  - Structure of input (observed) echos (FLASH [+ MT pulse])
% out - Structure of output (fitted) parameters (A, R1, R2, MT)
% pre - Structure of pre-computed parameters (B1t, B1r)
% opt - Structure of options
%
% OUTPUT
% ------
% llx - Data log-likelihood
% g   - Gradient of parameters
% H   - Hessian of parameters
% r   - Residuals

    set_num_threads(-1);
    
    if nargin < 5 || isempty(opt)
        opt = struct;
    end
    if nargin < 4 || isempty(pre)
        pre = struct;
    end
    
    % ---------------------------------------------------------------------
    % Read individual info
    dat = in.dat{n};                % File array
    v   = in.vol(n);                % Volume index
    TE  = in.TE(n);                 % Echo time
    TR  = in.TR(v);                 % Repetition time
    FA  = in.FA(v);                 % Flip angle (radian)
    mtpulse = in.mtpulse(v);        % Has MT pulse
    mat = in.mat(:,:,v);            % Observed-to-Common matrix
    dim = in.dim(:,v)';             % Observed dimensions
    lam = 1./in.var(v);             % Noise precision
    
    % ---------------------------------------------------------------------
    % Undersampling
    vs  = sqrt(sum(in.mat(1:3,1:3,v).^2)); % Voxel size
    if isfield(opt, 'subsample')
        skip = opt.subsample;
    else
        skip = [Inf Inf Inf];
    end
    skip(~isfinite(skip)) = vs(~isfinite(skip));
    skip  = round(skip ./ vs);
    skip  = max(skip, 1);
    Mskip = [diag(skip) (ones(3,1)-skip(:)); zeros(1,3) 1];
    mat   = mat * Mskip;
    dim   = ceil(dim ./ skip);
    
    % ---------------------------------------------------------------------
    % Parameter indices
    iA   = out.A.idx;
    iR1  = out.R1.idx;
    iR2  = out.R2.idx;
    iMT  = out.MT.idx;
    K    = isfinite(iA) + isfinite(iR1) + isfinite(iR2) + isfinite(iMT);
    
    % ---------------------------------------------------------------------
    % Load data
    
    % I am loading all the data at once to be fast. If it requires too much
    % memory, we could load it slice-wise inside the loop. But that means 
    % more i/o and would therefore be slower.
    
    % It is a shame that Matlab does not allow multi-threading with
    % shared memory, because the current estimates could be shared across
    % processes, really.
    
    lA  = pull(single(out.dat(:,:,:,iA)),  dim, out.mat\mat);
    lR1 = pull(single(out.dat(:,:,:,iR1)), dim, out.mat\mat);
    lR2 = pull(single(out.dat(:,:,:,iR2)), dim, out.mat\mat);
    if mtpulse
        lMT = pull(single(out.dat(:,:,:,iMT)), dim, out.mat\mat);
    end
    x   = single(dat(1:skip(1):end, ...
                     1:skip(2):end, ...
                     1:skip(3):end));
    lB1t = zeros([1 1 size(x,3)]);
    if isfield(pre, 'B1t')
        i_B1t = find(pre.B1t.vol == v);
        if isempty(i_B1t)
            i_B1t = find(pre.B1t.vol == Inf);
        end
        if ~isempty(i_B1t)
            lB1t = pull(single(pre.B1t.dat{i_B1t}()), dim, pre.B1t.mat(:,:,i_B1t)\mat);
            lB1t = log(max(lB1t/100, eps('single')));
        end
    end
    lB1r = zeros([1 1 size(x,3)]);
    if isfield(pre, 'B1r')
        i_B1r = find(pre.B1r.vol == v);
        if isempty(i_B1r)
            i_B1r = find(pre.B1r.vol == Inf);
        end
        if ~isempty(i_B1r)
            lB1r = pull(single(pre.B1r.dat{i_B1r}()), dim, pre.B1r.mat(:,:,i_B1r)\mat);
            lB1r = log(max(lB1r, eps('single')));
        end
    end
    
    % ---------------------------------------------------------------------
    % Allocate gradient/Hessian
    llx = 0;
    if nargout > 1
        g = zeros([dim K], 'single');
        if nargout > 2
            [ind, KK] = symIndices(K);
            H = zeros([dim KK], 'single');
        end
        if nargout > 3
            res = zeros(dim, 'single');
        end
    end
    
    % ---------------------------------------------------------------------
    % Loop over slices
    for z=1:size(x,3)
    
        % -----------------------------------------------------------------
        % Signal fit
        % * Proton density (x receive sensitivity)
        A  = exp(double(lA(:,:,z)) + double(lB1r(:,:,z)));
        % * T2 exponential decay: exp(-R2*TE)
        TER2 = exp(double(lR2(:,:,z)) + log(TE));   % decay value
        eR2  = exp(-TER2);                          % exponential term
        % * T1 exponential decay: exp(-R1*TR)
        TRR1  = exp(double(lR1(:,:,z)) + log(TR));  % decay value
        eR1   = exp(-TRR1);                         % exponential term
        % * True flip angle (x transmit sensitivity)
        sina = exp(log(FA) + double(lB1t(:,:,z)));
        cosa = 1 - 0.5 * sina.^2;
        % * MT ratio
        if mtpulse
            emt  = exp(-double(lMT(:,:,z)));
            d    = 1./(1 + emt);                 % MT-ratio
            md   = 1 - d;                        % 1 - MTratio
            dmd  = - emt./(1 + emt).^2;          % first derivative
            clear emt
        end
        % * Echo fit
        if mtpulse
            y = md.* (1 - eR1) ./ (1 - cosa .* md .* eR1);
        else
            y = (1 - eR1) ./ (1 - cosa .* eR1);
        end
        y = y .* A .* sina .* eR2;
        % * Oberved echo
        xz  = x(:,:,z);
        msk = isfinite(xz) & (xz > 0);
        % * Residuals
        r   = y - xz;
        clear xz
        
        % -----------------------------------------------------------------
        % Compute log-likelihood
        r(~isfinite(r)) = 0;
        llx = llx + sum(double(r(msk)).^2);

        % -----------------------------------------------------------------
        % Compute gradient and Hessian in observed space
        if nargout > 1
            % -------------------------------------------------------------
            % Proton density (A)
            if ~isnan(iA)
                g(:,:,z,iA) = y;
            end
            % -------------------------------------------------------------
            % T2 decay (R2)
            if ~isnan(iR2)
                g(:,:,z,iR2) = -TER2 .* y;
            end
            % -------------------------------------------------------------
            % T1 decay (R1)
            if ~isnan(iR1)
                if mtpulse
                    g(:,:,z,iR1) = md .* cosa ./ (1 - md .* cosa .* eR1) - 1 ./ (1 - eR1);
                else
                    g(:,:,z,iR1) = cosa ./ (1 - cosa .* eR1) - 1 ./ (1 - eR1);
                end
                g(:,:,z,iR1) = -TRR1 .* eR1 .* y .* double(g(:,:,z,iR1));
            end
            % -------------------------------------------------------------
            % MT ratio (MT)
            if ~isnan(iMT) && mtpulse
                g(:,:,z,iMT) = 1 ./ (md .* (1 - md .* cosa .* eR1));
                g(:,:,z,iMT) = dmd .* y .* double(g(:,:,z,iMT));
            end
        end
        if nargout > 2
            % -------------------------------------------------------------
            % The expectation of the Hessian is obtained from the gradient
            % (squared), up to the residuals.
            for k=1:K
                for kk=k:K
                    H(:,:,z,ind(k,kk)) = double(g(:,:,z,k)) .* double(g(:,:,z,kk));
                end
            end
        end
        if nargout > 1
            % -------------------------------------------------------------
            % Multiply gradient with residuals
            g(:,:,z,:) = bsxfun(@times, g(:,:,z,:), r);
        end
        if nargout > 3
            % -------------------------------------------------------------
            % Store residuals
            res(:,:,z) = r;
        end
        
    end
    
    % ---------------------------------------------------------------------
    % Noise precision + push to common space
    llx = - 0.5 * lam * llx;
    if nargout > 1
        g = g * lam;
        g(~isfinite(g)) = 0;
        g = bsxfun(@times, g, msk);
        g = push(g, dim, out.mat\mat, out.dim);
        if nargout > 2
            H = H * lam;
            H(~isfinite(H)) = 0;
            H = bsxfun(@times, H, msk);
            H = push(H, dim, out.mat\mat, out.dim);
        end
    end
    
    
    set_num_threads(1);
end

% -------------------------------------------------------------------------
% Check gradients for MPM nonlinear fit
% -------------------------------------------------------------------------
% 
% syms TE TR sina cosa 'positive'   % Fixed parameters
% syms lR1 lR2 lA ld 'real'         % Optimised parameters (log)
% syms X 'real'                     % Observation
% 
% % Precompute useful values
% A  = exp(lA);
% R1 = exp(lR1);
% R2 = exp(lR2);
% d  = 1 / (1 + exp(-ld));
% md = 1 - d;
% dd = - exp(-ld)./(1 + exp(-ld)).^2; % Derivative of (1-d) w.r.t. ld
% e1 = exp(-TR*R1);
% e2 = exp(-TE*R2);
% 
% % Signal fit and negative log-likelihood
% S = A * sina * md * (1 - e1) / (1 - md * cosa * e1) * e2;
% R = S - X;
% L = 0.5*R^2;
% 
% % Compute gradients automatically
% g = [diff(L, lA); ...
%      diff(L, lR2); ...
%      diff(L, lR1); ...
%      diff(L, ld)];
% H = [diff(L, lA,  lA) diff(L, lA,  lR2) diff(L, lA,  lR1) diff(L, lA,  ld); ...
%      diff(L, lR2, lA) diff(L, lR2, lR2) diff(L, lR2, lR1) diff(L, lR2, ld); ...
%      diff(L, lR1, lA) diff(L, lR1, lR2) diff(L, lR1, lR1) diff(L, lR1, ld); ...
%      diff(L, ld,  lA) diff(L, ld,  lR2) diff(L, ld,  lR1) diff(L, ld,  ld)];
% H = subs(H, X, S);
% 
% % Check that our gradients are correct
% gg = [S; ...
%       -TE * R2 * S; ...
%       -TR * R1 * e1 * S * (md * cosa / (1 - md * cosa .* e1) - 1 ./ (1 - e1)); ...
%       dd * S / (md * (1 - md * cosa * e1))];
% HH = gg * gg.';
% gg = gg * R;
% 
% simplify(g-gg, 100)
% simplify(H-HH, 100)