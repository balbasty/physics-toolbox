function llm = conditional(coils, sens, mean, prec, mask, senslog)
% Compute conditional log-likelihood
%
% FORMAT ll = b1m.ll.conditional(coils, sens, mean, prep, [mask], [senslog])
%
% coils   - (File)Array [Nx Ny Nz Nc] - Complex coil images
% sens    - (File)Array [Nx Ny Nz Nc] - Complex log-sensitivity profiles
% mean    - (File)Array [Nx Ny Nz]    - Complex mean image
% prec    -       Array [Nc Nc]       - Noise precision matrix
% mask    -       Array [Nx Ny]       - K-space sampling mask [full]
% senslog -                           - Log-sensitivities [false]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 6
    senslog = false;
    if nargin < 5
        mask = [];
    end
end

Nx = size(coils,1);
Ny = size(coils,2);
Nz = size(coils,3);
Nc = size(coils,4);
Nvox = Nx*Ny*Nz;

gpu_on = isa(prec, 'gpuArray');
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end

% -------------------------------------------------------------------------
% Compute log-likelihood (conditional)
llm = 0;
for z=1:Nz

    % ---------------------------------------------------------------------
    % Load one slice of the complete coil dataset
    xz = loadarray(coils(:,:,z,:), @single);
    xz = reshape(xz, [], Nc);

    % ---------------------------------------------------------------------
    % Compute map of missing data
    cz = utils.gmm.lib('obs2code', xz);
    code_list = unique(cz)';
    
    % ---------------------------------------------------------------------
    % Load one slice of the (previous) mean
    rz = loadarray(mean(:,:,z,:), @single);
    rz = reshape(rz, [], 1);

    % ---------------------------------------------------------------------
    % Load one slice of the complete sensitivity dataset + correct
    sz = loadarray(sens(:,:,z,:), @single);
    sz = reshape(sz, [], Nc);
    if senslog, sz = exp(sz); end
    rz = bsxfun(@times, rz, sz);

    % ---------------------------------------------------------------------
    % If incomplete sampling: push+pull coil-specific means
    prz = b1m.adjoint_forward(reshape(rz, [Nx Ny 1 Nc]), mask);
    prz = reshape(prz, [], Nc);

    % ---------------------------------------------------------------------
    % Missing data
    for code=code_list
        sub = cz == code;
        bin  = utils.gmm.lib('code2bin', code, Nc);
        if ~any(bin)
            continue
        end
        A = utils.invPD(prec);
        A = A(bin,bin);
        A = utils.invPD(A);
        
        llm = llm - sum(double(real(dot(rz(sub,bin),(prz(sub,bin)-2*xz(sub,bin))*A,2))));
    end
    

end
llm = 0.5 * Nvox * llm;