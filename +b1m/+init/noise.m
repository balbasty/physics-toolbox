function [C,A] = noise(coils, gpu)
% Estimate the noise variance of each coil by fitting a Rice or Gaussian
% mixture model to their magnitude image.
%
% FORMAT [cov,prec] = b1m.init.noise(coils, (gpu))
%
% coils - (File)Array [Nx Ny Nz Nc] - Complex coil images
% gpu                               - If 'gpu', compute on GPU
%
% cov   -       Array [Nc Nc]       - Noise covariance matrix
% prec  -       Array [Nc Nc]       - Noise precision matrix
%
% >> C = sum((s*r-x)*(s*r-x)')/(2J)     [where ' = conjugate transpose]
%
% Nc = number of coils
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 4
    gpu = '';
end
gpu_on = strcmpi(gpu, 'gpu');

% -------------------------------------------------------------------------
% Allocate output array
C = zeros(size(coils,4), 'single');
if gpu_on
    C = gpuArray(C);
end
if gpu_on, loadarray = @utils.loadarray_gpu;
else,      loadarray = @utils.loadarray_cpu; end

Nc = size(coils,4);
ndigits = ceil(log10(Nc))+1;

fprintf(['InitCov: %' num2str(ndigits) 's'], '');
for n=1:Nc
    
    fprintf([repmat('\b', [1 ndigits]) '%' num2str(ndigits) 'd'], n);
    
    % ---------------------------------------------------------------------
    % Load one coil image
    xn = loadarray(abs(coils(:,:,:,n)), @single);
    xn = reshape(xn, [], 1);
    
    % ---------------------------------------------------------------------
    % Compute histogram
    xmax = max(xn,[],'omitnan');
    xmin = min(xn,[],'omitnan');
    cn = linspace(xmin,xmax,129);
    cn = (cn(2:end)+cn(1:end-1))/2;
    bwn = (xmax - xmin)/128;
    [xn,wn] = spm_imbasics('hist', double(xn), cn(:), 'KeepZero', false);
    clear c1
    
    % ---------------------------------------------------------------------
    % Fit Rician mixture
    [PI,NU,SIG] = spm_rice_mixture(double(wn), double(xn), 2);
    if NU(1) == 0 && NU(2) == 0
        % If failed, fit Gaussian mixture
        [~,MU,A,PI] = utils.gmm.fit(double(xn), 2, double(wn),...
            'GaussPrior', {[],10,[],10},'BinWidth',bwn);
        if MU(1) <= MU(2)
            C(n,n) = 1./A(1);
        else
            C(n,n) = 1./A(2);
        end
    else
        if NU(1) <= NU(2)
            C(n,n) = SIG(1)^2;
        else
            C(n,n) = SIG(2)^2;
        end
    end
    
end
fprintf('\n');

A = utils.matcomp('inv', C);