function [s,llm,llp,ok] = multicoil_sensitivity(list_n, rho, x, s, A, prm, vs, llp, centre_fields)
% Maximum a posteriori sensitivity profiles given a set of observed coil 
% images, a mean image and a noise precision (= inverse covariance) 
% matrix.
%
% FORMAT s = multicoil_sensitivity(n, rho, x, s, A, prm, (vs))
%
% n   -                               - Index of the coil to update
% rho - (File)Array [Nx Ny Nz  1 (2)] - Complex mean image
% x   - (File)Array [Nx Ny Nz Nc (2)] - Complex coil images
% s   - (File)Array [Nx Ny Nz Nc (2)] - Complex log-sensitivity profiles
% A   -       Array [Nc Nc]           - Noise precision matrix
% prm -       Array [1 3] or [Nc 3]   - Regularisation (/ coil) [a m b]
% vs  -       Array [1 3]             - Voxel size [1 1 1]
%
% The optimum is found numerically using complex Gauss-Newton optimisation.
% The inverse problem is real and is solved by full multigrid.
%
% Nc = number of coils
% Images can either be complex or have two real components that are then 
% assumed to be the real and imaginary parts.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Neumann boundary condition (null derivative)
spm_field('boundary', 1); 

% -------------------------------------------------------------------------
% Default values
if isempty(list_n)
    list_n = 1:size(x,4); % Coils to process
end
list_n = list_n(:).'; % Ensure row-vector
if nargin < 7
    vs = [1 1 1]; % Voxel size
end
if nargin < 8
    llp = NaN;  % Previous log-likelihood (prior part)
end
if nargin < 9
    centre_fields = false; % Ensure that bias fields sum to zero
end
if size(prm, 1) == 1
    prm = repmat(prm, [size(x,4) 1]); % Use same parameters for all coils
end
% Compute coil-wise regularisation modulation
% We assume that prm(n) = prm_global * alpha(n), such that sum(alpha) = 1
alpha = sum(prm,2);
alpha = alpha/sum(alpha);
alpha = reshape(alpha, 1, []);
common_prm = prm(1,:) / alpha(1);
=
N = size(x,4);

% -------------------------------------------------------------------------
% For each coil
fprintf('Update Sensitivity:');
for n=list_n
    
    fprintf(' %2d', n);
    reg = [vs prm(n,:)];
    
    % ---------------------------------------------------------------------
    % Compute log-likelihood (prior)
    if isnan(llp)
        llp = multicoil_ll_prior(s, prm, vs);
    end
    % ---------------------------------------------------------------------
    % Prepare weights: beta(n,m) = (n == m) - alpha(n)
    beta    = repmat(-alpha(n), [1 N]);
    beta(n) = 1 + beta(n);
    
    % ---------------------------------------------------------------------
    % Allocate conjugate gradient and Hessian
    g = zeros(size(x,1), size(x,2), size(x,3), 2, 'single');
    H = zeros(size(x,1), size(x,2), size(x,3), 1, 'single');
    llm = 0;
    
    % ---------------------------------------------------------------------
    % Compute gradient slice-wise to save memory
    for z=1:size(rho, 3) 

        % -----------------------------------------------------------------
        % Load one slice of the complete coil dataset
        if size(x, 5) == 2
            % Two real components
            x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
            x1 = x1(:,:,1) + 1i*x1(:,:,2);
        else
            % One complex volume
            x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
        end
        if isa(A, 'gpuArray')
            x1 = gpuArray(x1);
        end

        % -----------------------------------------------------------------
        % Load one slice of the mean
        if size(rho, 5) == 2
            % Two real components
            rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
            rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
        else
            % One complex volume
            rho1 = reshape(single(rho(:,:,z,:,:)), [], 1);
        end
        if isa(A, 'gpuArray')
            rho1 = gpuArray(rho1);
        end

        % -----------------------------------------------------------------
        % Load one slice of the complete sensitivity dataset + correct
        if size(s, 5) == 2
            % Two real components
            s1 = reshape(double(s(:,:,z,:,:)), [], size(x,4), 2);
            s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
        else
            % One complex volume
            s1 = reshape(single(exp(double(s(:,:,z,:,:)))), [], size(x,4));
        end
        if isa(A, 'gpuArray')
            s1 = gpuArray(s1);
        end
        rho1 = bsxfun(@times, rho1, s1);
        clear s1
        
        % -----------------------------------------------------------------
        % Compute gradient
        
        tmp = (rho1 - x1) * A;
        
        llm = llm - 0.5 * sum(double(real(dot(tmp, rho1 - x1, 2))));
        
        if centre_fields
            rho1 = bsxfun(@times, rho1, beta);
            tmp  = dot(tmp, rho1, 2);
        else
            tmp = rho1(:,n) .* conj(tmp(:,n));
        end
        
        g(:,:,z,1) = g(:,:,z,1) + reshape( real(tmp), size(g,1), size(g,2));
        g(:,:,z,2) = g(:,:,z,2) + reshape(-imag(tmp), size(g,1), size(g,2));
        
        if centre_fields
            tmp = real(dot(rho1, rho1 * A, 2));
        else
            tmp = A(n,n) * real(conj(rho1(:,n)) .* rho1(:,n));
        end
        H(:,:,z) = H(:,:,z) + reshape(tmp, size(H,1), size(H,2));
        
        clear tmp
    end
    clear x1 rho1 g1 H1
    
    % ---------------------------------------------------------------------
    % Gather gradient & Hessian (if on GPU)
    g = gather(g);
    H = gather(H);
    
    % ---------------------------------------------------------------------
    % Load previous bias field
    if size(s, 5) == 2
        % Two real components
        s0 = single(s(:,:,:,n,:));
    else
        % One complex volume
        s0 = single(s(:,:,:,n));
        s0 = cat(5, real(s0), imag(s0));
    end
    
    % ---------------------------------------------------------------------
    % Prior term
    llpr = spm_field('vel2mom', s0(:,:,:,:,1), reg);
    g(:,:,:,1) = g(:,:,:,1) + llpr;
    llpi = spm_field('vel2mom', s0(:,:,:,:,2), reg);
    g(:,:,:,2) = g(:,:,:,2) + llpi;
    clear s0

    % ---------------------------------------------------------------------
    % Gauss-Newton
    regH = reg;
    if centre_fields
        regH(4:end) = regH(4:end) * beta(n);
    end
    ds1 = spm_field(H, g(:,:,:,1), [regH 2 2]);
    ds2 = spm_field(H, g(:,:,:,2), [regH 2 2]);
    clear g H
    ds = cat(5, ds1, ds2);
    clear ds1 ds2
    
    % ---------------------------------------------------------------------
    % Parts for log-likelihood (prior)
    Lds = cat(4, spm_field('vel2mom', ds(:,:,:,1), [vs common_prm]), ...
                 spm_field('vel2mom', ds(:,:,:,2), [vs common_prm]));
    if centre_fields
        sums = 0;
        sumb = sum(alpha .* beta(:)'.^2);
        for m=1:N
            if size(s, 5) == 2
                % Two real components
                s1 = single(s(:,:,:,m,:));
            else
                % One complex volume
                s1 = single(s(:,:,:,m,:));
                s1 = cat(5, real(s1), imag(s1));
            end
            sums = sums + alpha(m) * beta(m) * s1;
            clear s1
        end
        % part1 = (sum alpha*beta) * (ds)'L(ds)
        % part2 = -2 * (sum alpha*beta*s)'L(ds)
        part1r = sumb * double(reshape(ds(:,:,:,1), 1, [])) * double(reshape(Lds(:,:,:,1), [], 1));
        part1i = sumb * double(reshape(ds(:,:,:,2), 1, [])) * double(reshape(Lds(:,:,:,2), [], 1));
        part2r = -2 * double(reshape(Lds(:,:,:,1), 1, [])) * double(reshape(sums(:,:,:,1), [], 1));
        part2i = -2 * double(reshape(Lds(:,:,:,2), 1, [])) * double(reshape(sums(:,:,:,2), [], 1));
        clear Lds sums
    else
        partr = alpha(n) * double(reshape(ds(:,:,:,1), 1, [])) * double(reshape(Lds(:,:,:,1), [], 1));
        parti = alpha(n) * double(reshape(ds(:,:,:,2), 1, [])) * double(reshape(Lds(:,:,:,2), [], 1));
        clear Lds
    end
    
    beta = reshape(beta, [1 1 1 N]);
        
    % ---------------------------------------------------------------------
    % Line-Search
    llm0   = llm;
    llp0   = llp;
    ok     = false;
    armijo = 1;
    for ls=1:6
        
        % -----------------------------------------------------------------
        % Compute log-likelihood (prior)
        if centre_fields
            llpr = -0.5 * (armijo^2 * part1r + armijo * part2r);
            llpi = -0.5 * (armijo^2 * part1i + armijo * part2i);
        else
            llpr = -0.5 * armijo^2 * partr;
            llpi = -0.5 * armijo^2 * parti;
        end
        llp  = llp0 + llpr + llpi;
        
        % -----------------------------------------------------------------
        % Compute log-likelihood (conditional)
        llm = 0;
        for z=1:size(rho, 3) 

            % -------------------------------------------------------------
            % Load one slice of the complete coil dataset
            if size(x, 5) == 2
                % Two real components
                x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
                x1 = x1(:,:,1) + 1i*x1(:,:,2);
            else
                % One complex volume
                x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
            end
            if isa(A, 'gpuArray')
                x1 = gpuArray(x1);
            end

            % -------------------------------------------------------------
            % Load one slice of the mean
            if size(rho, 5) == 2
                % Two real components
                rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
                rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
            else
                % One complex volume
                rho1 = reshape(single(rho(:,:,z,:,:)), [], 1);
            end
            if isa(A, 'gpuArray')
                rho1 = gpuArray(rho1);
            end

            % -------------------------------------------------------------
            % Load one slice of the complete sensitivity dataset + correct
            if size(s, 5) == 2
                % Two real components
                s1 = single(s(:,:,z,:,:));
                if centre_fields
                    s1(:,:,:,:,:) = s1(:,:,:,:,:) ...
                        - armijo * bsxfun(@times, beta, ds(:,:,z,:,:));
                else
                    s1(:,:,:,n,:) = s1(:,:,:,n,:) - armijo * ds(:,:,z,:,:);
                end
                s1 = reshape(double(s1), [], size(x,4), 2);
                s1 = single(exp(s1(:,:,1) + 1i*s1(:,:,2)));
            else
                % One complex volume
                s1 = single(s(:,:,z,:,:));
                if centre_fields
                    s1(:,:,:,:,:) = s1(:,:,:,:,:) ...
                        -  armijo * bsxfun(@times, beta, ds(:,:,z,:,1) + 1i * ds(:,:,z,:,2));
                else
                    s1(:,:,:,n,:) = s1(:,:,:,n,:) - armijo * (ds(:,:,z,:,1) + 1i * ds(:,:,z,:,2));
                end
                s1 = reshape(single(exp(double(s1))), [], size(x,4));
            end
            if isa(A, 'gpuArray')
                s1 = gpuArray(s1);
            end
            rho1 = bsxfun(@times, rho1, s1);
            clear s1

            % -------------------------------------------------------------
            % Compute gradient

            llm = llm - 0.5 * sum(double(real(dot((rho1 - x1) * A, rho1 - x1, 2))));

        end
        clear x1 rho1 g1 H1
        
        % -----------------------------------------------------------------
        % Check progress
        if (llm+llp) > (llm0+llp0)
            ok = true;
            break
        else
            armijo = armijo/2;
        end
        
    end
    
    % ---------------------------------------------------------------------
    % Save
    if ok
        fprintf(' :D (%d)', ls);
        if size(s,5) == 1
            ds = ds(:,:,:,1) + 1i * ds(:,:,:,2);
        end
        if centre_fields
            for m=1:size(x,4)
                s(:,:,:,m,:) = s(:,:,:,m,:) - beta(m) * armijo * ds;
            end
        else
            s(:,:,:,n,:) = s(:,:,:,n,:) - armijo * ds;
        end
    else
        fprintf(' :(');
        llm = llm0;
        llp = llp0;
    end
        
end
fprintf('\n');