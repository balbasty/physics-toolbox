function s = multicoil_init_phase(rho, x, s)

% -------------------------------------------------------------------------
% Allocate output volume
if nargin < 3
    s = [];
end
if isempty(s)
    s = zeros(size(x,1), size(x,2), size(x,3), size(x,4), 2, 'single');
end

% -------------------------------------------------------------------------
% Process slice-wise to save memory
sumsin = zeros(1,size(x,4));
sumcos = zeros(1,size(x,4));
for z=1:size(rho, 3) 

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    if size(x, 5) == 2
        % Two real components
        x1 = reshape(single(x(:,:,z,:,:)), [], size(x,4), 2);
        x1 = x1(:,:,1) + 1i*x1(:,:,2);
    else
        % One complex volume
        x1 = reshape(single(x(:,:,z,:)), [], size(x,4));
    end
%     if isa(A, 'gpuArray')
%         x1 = gpuArray(x1);
%     end
    x1 = angle(x1);

    
    % ---------------------------------------------------------------------
    % Load one slice of the mean
    if size(rho, 5) == 2
        % Two real components
        rho1 = reshape(single(rho(:,:,z,:,:)), [], 1, 2);
        rho1 = rho1(:,:,1) + 1i*rho1(:,:,2);
    else
        % One complex volume
        rho1 = reshape(single(rho(:,:,z,:)), [], 1);
    end
%     if isa(A, 'gpuArray')
%         rho1 = gpuArray(rho1);
%     end
    rho1 = angle(rho1);

    % ---------------------------------------------------------------------
    % Accumulate sin(diff) and cos(diff)
    s1 = bsxfun(@minus, x1, rho1);
    sumsin = sumsin + sum(sin(s1),1);
    sumcos = sumcos + sum(cos(s1),1);

end

% -------------------------------------------------------------------------
% Compute mean shift (Von Mises maximum likelihood)
shift = atan2(sumsin,sumcos);

% ---------------------------------------------------------------------
% Write
for n=1:size(s,4)
    if size(s, 5) == 2
        s(:,:,:,n,2) = shift(n);
    else
        s(:,:,:,n) = s(:,:,:,n) - 1i*imag(s(:,:,:,n)) + 1i*shift(n);
    end
end